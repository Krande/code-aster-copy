! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine lcotan(opt, etatd, etatf, fami, &
                  kpg, ksp, rela_comp, mod, imat, &
                  nmat, materd, materf, epsd, deps, &
                  sigd, sigf, nvi, vind, vinf, &
                  drdy, vp, vecp, theta, dt, &
                  devg, devgii, timed, timef, compor, &
                  nbcomm, cpmono, pgl, nfs, nsg, &
                  toutms, hsr, nr, itmax, toler, &
                  typma, dsde, codret)
! aslint: disable=W1504
    implicit none
! ======================================================================
!
!     CALCUL DE L'OPERATEUR TANGENT = DS/DE(T+DT) OU DS/DE(T)
!     CONVENTION :
!                 SUFFIXE D : DEBUT DU PAS DE TEMPS
!                 SUFFIXE F : FIN DU PAS DE TEMPS
!     ==================================================================
!     ARGUMENTS
!
!     IN FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!        KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!        IMAT    ADRESSE DU MATERIAU CODE
!        COMPOR  COMPORTEMENT DE L ELEMENT
!                COMP(1) = RELATION DE COMPORTEMENT (ROUSSELIER.)
!                COMP(2) = NB DE VARIABLES INTERNES
!                COMP(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!        TIMED   INSTANT T
!        TIMEF   INSTANT T+DT
!        TEMPD   TEMPERATURE A T           POUR LA THM
!        TEMPF   TEMPERATURE A T+DT        POUR LA THM
!        TREF    TEMPERATURE DE REFERENCE  POUR LA THM
!        CES PARAMETRES DE TEMPERATURE NE SONT PAS PRIS EN COMPTE EN
!        MECANIQUE PURE (ON UTILISE LES VARIABLES DE COMMANDES)
!
!        EPSDT   DEFORMATION TOTALE A T
!        DEPST   INCREMENT DE DEFORMATION TOTALE
!        SIGD    CONTRAINTE A T
!        VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!        OPT     OPTION DE CALCUL
!                        'RIGI_MECA_TANG'> DSDE(T)
!                        'FULL_MECA'     > DSDE(T+DT), SIGF, VINF
!                        'RAPH_MECA'     > SIGF, VINF
!     OUT
!        SIGF    CONTRAINTE A T+DT
!        VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!        DSDE    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!        CODRET  CODE RETOUR =0 OK, =1 => REDECOUPAGE DU PAS DE TEMPS
!
#include "asterfort/lchbvp.h"
#include "asterfort/lcjela.h"
#include "asterfort/lcjpla.h"
#include "asterfort/lcjplc.h"
#include "asterfort/Behaviour_type.h"
    integer(kind=8) :: nmat, nsg, nfs, nbcomm(nmat, 3)
!
    character(len=*) :: fami
    character(len=7) :: etatd, etatf
    character(len=8) :: mod, typma
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=16) :: opt
    character(len=24) :: cpmono(5*nmat+1)
!
    integer(kind=8) :: imat, ndt, ndi, nr, nvi, itmax, kpg, ksp, codret, k, l
    real(kind=8) :: toler, materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: vind(*), vinf(*), timed, timef, epsd(9), deps(9), sigd(6)
    real(kind=8) :: sigf(6)
    real(kind=8) :: theta, dt, devg(6), devgii, vp(3), vecp(3, 3), dsde(6, *)
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg), drdy(*)
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
!     OPTIONS 'FULL_MECA' ET 'RIGI_MECA_TANG' = CALCUL DE DSDE
!     ----------------------------------------------------------------
!     EVALUATION DU JACOBIEN DSDE A (T+DT) POUR 'FULL_MECA'
!     ET CALCUL ELASTIQUE    ET   A (T)    POUR 'RIGI_MECA_TANG'
!     ----------------------------------------------------------------
!
    codret = 0
!
!     MATRICE TANGENTE DE PREDICTION
!
    if (opt(1:9) .eq. 'RIGI_MECA') then
!
        if (rela_comp .eq. 'LAIGLE') then
!
            call lcjela(rela_comp, mod, nmat, materd, vind, &
                        dsde)
!
        else if ((etatd .eq. 'PLASTIC') .and. (rela_comp .eq. 'MONOCRISTAL')) then
!
            call lcjplc(rela_comp, mod, nmat, &
                        materf, timed, timef, compor, nbcomm, &
                        cpmono, pgl, nfs, nsg, toutms, &
                        hsr, nr, nvi, epsd, deps, &
                        itmax, toler, sigd, vind, sigd, &
                        vind, dsde, drdy, opt, codret)
            if (codret .ne. 0) goto 999
!
        else if ((etatd .eq. 'PLASTIC') .and. (typma .eq. 'VITESSE ')) then
            if ((rela_comp .eq. 'HOEK_BROWN') .or. (rela_comp .eq. 'HOEK_BROWN_EFF')) then
! ---              HOEK-BROWN : CALCUL DES VALEURS ET VECTEURS PROPRES
! ---                           DU DEVIATEUR ELASTIQUE
                call lchbvp(sigd, vp, vecp)
            end if
            call lcjpla(fami, kpg, ksp, rela_comp, mod, &
                        nr, imat, nmat, materd, nvi, &
                        deps, sigd, vind, dsde, sigd, &
                        vind, vp, vecp, theta, dt, &
                        devg, devgii, codret)
            if (codret .ne. 0) goto 999
!
        else
!
!           CAS GENERAL : ELASTICITE LINEAIRE ISOTROPE OU ANISOTROPE
            call lcjela(rela_comp, mod, nmat, materd, vind, &
                        dsde)
!
        end if
!
!
    else if (opt(1:9) .eq. 'FULL_MECA') then
!
        if (etatf .eq. 'ELASTIC') then
            call lcjela(rela_comp, mod, nmat, materf, vinf, dsde)
!
        else if (etatf .eq. 'PLASTIC') then
!   --->    ELASTOPLASTICITE ==>  TYPMA = 'VITESSE '
!   --->    VISCOPLASTICITE  ==>  TYPMA = 'COHERENT '
            if (typma .eq. 'COHERENT') then
                call lcjplc(rela_comp, mod, nmat, &
                            materf, timed, timef, compor, nbcomm, &
                            cpmono, pgl, nfs, nsg, toutms, &
                            hsr, nr, nvi, epsd, deps, &
                            itmax, toler, sigf, vinf, sigd, &
                            vind, dsde, drdy, opt, codret)
                if (codret .ne. 0) goto 999
            else if (typma .eq. 'VITESSE ') then
                call lcjpla(fami, kpg, ksp, rela_comp, mod, &
                            nr, imat, nmat, materd, nvi, &
                            deps, sigf, vinf, dsde, sigd, &
                            vind, vp, vecp, theta, dt, &
                            devg, devgii, codret)
                if (codret .ne. 0) goto 999
            end if
        end if
!
    end if
!
! -   MODIFICATION EN CONTRAINTE PLANES POUR TENIR COMPTE DE
!     SIG3=0 ET DE LA CONSERVATION DE L'ENERGIE
    if (mod(1:6) .eq. 'C_PLAN') then
        do k = 1, ndt
            if (k .eq. 3) cycle
            do l = 1, ndt
                if (l .eq. 3) cycle
                dsde(k, l) = dsde(k, l)-1.d0/dsde(3, 3)*dsde(k, 3)*dsde(3, l)
            end do
        end do
    end if
!
999 continue
end subroutine
