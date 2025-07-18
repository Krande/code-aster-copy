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
! aslint: disable=W1306,W1504
!
subroutine lcplnl(BEHinteg, &
                  fami, kpg, ksp, rela_comp, toler, &
                  itmax, mod, imat, nmat, materd, &
                  materf, nr, nvi, timed, timef, &
                  deps, epsd, sigd, vind, compor, &
                  nbcomm, cpmono, pgl, nfs, nsg, &
                  toutms, hsr, sigf, vinf, icomp, &
                  codret, drdy, carcri)
!
    use Behaviour_type
!
    implicit none
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
!
!     INTEGRATION ELASTO-PLASTIQUE ET VISCO-PLASTICITE
!           SUR DT DE Y = ( SIG , VIN )
!     LE SYSTEME  A RESOUDRE EN DY ETANT NON  LINEAIRE
!
!     ON RESOUD DONC                  R(DY) = 0
!     PAR UNE METHODE DE NEWTON       DRDY(DYI) DDYI = - R(DYI)
!                                     DYI+1 = DYI + DDYI  (DYO DEBUT)
!     ET ON REACTUALISE               YF = YD + DY
!
!     ATTENTION :     ON REACTUALISE ICI DEPS DE FACON A CE QUE
!                     DEPS(3) = DY(NR) EN C_PLAN
!
!
!     IN  FAMI   :  FAMILLE DE POINT DE GAUSS
!         KPG    :  NUMERO DU POINT DE GAUSS
!         KSP    :  NUMERO DU SOUS-POINT DE GAUSS
!         TOLER  :  TOLERANCE DE CONVERGENCE LOCALE
!         ITMAX  :  NOMBRE MAXI D'ITERATIONS LOCALES
!         MOD    :  TYPE DE MODELISATION
!         IMAT   :  ADRESSE DU MATERIAU CODE
!         NMAT   :  DIMENSION MATER
!         MATERD :  COEFFICIENTS MATERIAU A T
!         MATERF :  COEFFICIENTS MATERIAU A T+DT
!         NR     :  NB EQUATION DU SYSTEME R(DY)
!         NVI    :  NB VARIABLES INTERNES
!         TIMED  :  INSTANT  T
!         TIMEF  :  INSTANT T+DT
!         EPSD   :  DEFORMATION A T
!         SIGD   :  CONTRAINTE A T
!         VIND   :  VARIABLES INTERNES A T
!         COMP   :  COMPOR - LOI ET TYPE DE DEFORMATION
!         NBCOMM :  INCIDES DES COEF MATERIAU monocristal
!         CPMONO :  NOM DES COMPORTEMENTS monocristal
!         PGL    :  MATRICE DE PASSAGE
!         TOUTMS :  TENSEURS D'ORIENTATION monocristal
!         HSR    :  MATRICE D'INTERACTION monocristal
!         ICOMP  :  COMPTEUR POUR LE REDECOUPAGE DU PAS DE TEMPS
!     VAR DEPS   :  INCREMENT DE DEFORMATION
!     OUT SIGF   :  CONTRAINTE A T+DT
!         VINF   :  VARIABLES INTERNES A T+DT
!         CODRET :  CONTROLE DU REDECOUPAGE DU PAS DE TEMPS
!         DRDY   :  JACOBIEN
!         CRIT   :  CRITERES DE CONVERGENCE LOCAUX
! VARIABLES LOCALES
!         R      :  VECTEUR RESIDU
!         DY     :  INCREMENT DES VARIABLES = ( DSIG  DVIN  (DEPS3)  )
!         DDY    :  CORRECTION SUR L'INCREMENT DES VARIABLES
!                                           = ( DDSIG DDVIN (DDEPS3) )
!         YD     :  VARIABLES A T   = ( SIGD  VARD  )
!         YF     :  VARIABLES A T+DT= ( SIGF  VARF  )
!         TYPESS :  TYPE DE SOLUTION D ESSAI POUR NEWTON
!         ESSAI  :  VALEUR  SOLUTION D ESSAI POUR NEWTON
!         INTG   :  COMPTEUR DU NOMBRE DE TENTATIVES D'INTEGRATIONS
!
#include "asterf_types.h"
#include "asterfort/lcafyd.h"
#include "asterfort/lccaga.h"
#include "asterfort/lcconv.h"
#include "asterfort/lcinit.h"
#include "asterfort/lcjacb.h"
#include "asterfort/lcjacp.h"
#include "asterfort/lcplnf.h"
#include "asterfort/lcreli.h"
#include "asterfort/lcresi.h"
#include "asterfort/mgauss.h"
#include "asterfort/r8inir.h"
#include "asterfort/utlcal.h"
#include "asterfort/Behaviour_type.h"
    integer(kind=8) :: imat, nmat, icomp
    integer(kind=8) :: typess, itmax, iret, kpg, ksp
    integer(kind=8) :: nr, ndt, ndi, nvi, iter
!
    real(kind=8) :: toler, essai, rbid
    real(kind=8) :: epsd(*), deps(*)
    real(kind=8) :: sigd(6), sigf(6)
    real(kind=8) :: vind(*), vinf(*)
!     DIMENSIONNEMENT DYNAMIQUE (MERCI F90)
    real(kind=8) :: r(nr), drdy(nr, nr), rini(nr)
    real(kind=8) :: drdy1(nr, nr)
    real(kind=8) :: ddy(ndt+nvi), dy(ndt+nvi), yd(ndt+nvi), yf(ndt+nvi)
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), dt
    real(kind=8) :: timed, timef, drdyb(nr, nr)
    aster_logical :: lreli
!
    character(len=8) :: mod
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=*) :: fami
!
    common/tdim/ndt, ndi
!
    integer(kind=8) :: intg, codret
!
    integer(kind=8) :: nbcomm(nmat, 3), verjac, nfs, nsg
    real(kind=8) :: pgl(3, 3), epstr(6)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
    character(len=4) :: cargau
    character(len=16) :: algo
    character(len=24) :: cpmono(5*nmat+1)
!
    integer(kind=8) :: nr1
!
!     ACTIVATION OU PAS DE LA RECHERCHE LINEAIRE
    lreli = .false.
    call utlcal('VALE_NOM', algo, carcri(6))
    if (algo .eq. 'NEWTON_RELI') lreli = .true.
!
!     VERIFICATION DE LA MATRICE JACOBIENNE
!     VERJAC=0 : PAS DE VERIFICATION
!     VERJAC=1 : CONSTRUCTION DE LA JACOBIENNE PAR PERTURBATION (LCJACP)
!                COMPARAISON A LA MATRICE JACOBIENNE ISSU DE LCJACB
!     VERJAC=2 : UTILISATION DE LA JACOBIENNE PAR PERTURBATION (LCJACP)
!                COMME MATRICE JACOBIENNE A LA PLACE DE LCJACB
!
    verjac = 0
!
    if (algo .eq. 'NEWTON_PERT') then
        verjac = 2
    end if
!
    essai = 1.d-5
    dt = timef-timed
!
    call r8inir(ndt+nvi, 0.d0, yd, 1)
!
    codret = 0
!
    nr1 = nr
    iret = 0
!
!
!     CHOIX DES VALEURS DE VIND A AFFECTER A YD
    call lcafyd(compor, materd, materf, nbcomm, cpmono, &
                nmat, mod, nvi, vind, &
                sigd, nr1, yd)
!
!     CHOIX DES PARAMETRES DE LANCEMENT DE MGAUSS
    call lccaga(rela_comp, cargau)
!
    if (mod(1:6) .eq. 'C_PLAN') yd(nr) = epsd(3)
!
!     RESOLUTION ITERATIVE PAR NEWTON DE R(DY) = 0
!     SOIT  DRDY(DYI) DDYI = -R(DYI)  ET DYI+1 = DYI + DDYI
!                                                         -
! -   INITIALISATION DU TYPE DE SOLUTION D ESSAI (-1)
    typess = -1
    intg = 0
!
2   continue
!
    call r8inir(nr, 0.d0, r, 1)
    call r8inir(ndt+nvi, 0.d0, ddy, 1)
    call r8inir(ndt+nvi, 0.d0, dy, 1)
    call r8inir(ndt+nvi, 0.d0, yf, 1)
!
!     CALCUL DE LA SOLUTION D ESSAI INITIALE DU SYSTEME NL EN DY
    call lcinit(fami, kpg, ksp, rela_comp, typess, &
                essai, mod, nmat, materf, &
                timed, timef, nr1, nvi, &
                yd, epsd, deps, dy, compor, &
                nbcomm, cpmono, pgl, nfs, nsg, &
                toutms, vind, sigd, sigf, epstr, &
                iret)
!
    if (iret .ne. 0) then
        goto 3
    end if
!
    iter = 0
!
1   continue
!
!     ITERATIONS DE NEWTON
    iter = iter+1
!
!     PAR SOUCIS DE PERFORMANCES, ON NE REFAIT PAS DES OPERATIONS
!     QUI ONT DEJA ETE FAITE A L'ITERATION PRECEDENTE DANS LE CAS
!     DE LA RECHERCHE LINEAIRE
    if (.not. lreli .or. iter .eq. 1) then
!        INCREMENTATION DE  YF = YD + DY
        yf(1:nr) = yd(1:nr)+dy(1:nr)
!
!        CALCUL DES TERMES DU SYSTEME A T+DT = -R(DY)
        call lcresi(fami, kpg, ksp, rela_comp, mod, &
                    imat, nmat, materd, materf, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, nr, nvi, vind, &
                    vinf, itmax, toler, timed, timef, &
                    yd, yf, deps, epsd, dy, &
                    r, iret, carcri)
!
        if (iret .ne. 0) then
            goto 3
        end if
!
    end if
!     SAUVEGARDE DE R(DY0) POUR TEST DE CONVERGENCE
    if (iter .eq. 1) then
        rini = r
    end if
!
!
    if (verjac .ne. 2) then
!         CALCUL DU JACOBIEN DU SYSTEME A T+DT = DRDY(DY)
        call lcjacb(fami, kpg, ksp, rela_comp, mod, &
                    nmat, materf, timed, timef, &
                    yf, deps, itmax, toler, nbcomm, &
                    cpmono, pgl, nfs, nsg, toutms, &
                    hsr, nr, nvi, vind, &
                    vinf, epsd, yd, dy, &
                    carcri, &
                    drdy, iret)
        if (iret .ne. 0) then
            goto 3
        end if
    end if
!
    if (verjac .ge. 1) then
        call lcjacp(fami, kpg, ksp, rela_comp, toler, &
                    itmax, mod, imat, nmat, materd, &
                    materf, nr, nvi, timed, timef, &
                    deps, epsd, vind, vinf, yd, &
                    nbcomm, cpmono, pgl, nfs, &
                    nsg, toutms, hsr, dy, r, &
                    drdy, verjac, drdyb, iret, carcri)
        if (iret .ne. 0) goto 3
    end if
!
!     RESOLUTION DU SYSTEME LINEAIRE DRDY(DY).DDY = -R(DY)
    drdy1 = drdy
    ddy(1:nr) = r(1:nr)
    call mgauss(cargau, drdy1, ddy, nr, nr1, &
                1, rbid, iret)
!
    if (iret .ne. 0) then
        goto 3
    end if
!
!     ACTUALISATION DE LA SOLUTION
    if (.not. lreli) then
!        REACTUALISATION DE DY = DY + DDY
        dy(1:nr) = ddy(1:nr)+dy(1:nr)
    else if (lreli) then
!        RECHERCHE LINEAIRE : RENVOIE DY, YF ET R RE-ACTUALISES
        call lcreli(fami, kpg, ksp, rela_comp, mod, &
                    imat, nmat, materd, materf, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, nr, nvi, vind, &
                    vinf, itmax, toler, timed, timef, &
                    yd, yf, deps, epsd, dy, &
                    r, ddy, iret, carcri)
        if (iret .ne. 0) goto 3
    end if
    if (mod(1:6) .eq. 'C_PLAN') deps(3) = dy(nr)
!
!     VERIFICATION DE LA CONVERGENCE EN DY  ET RE-INTEGRATION ?
    call lcconv(rela_comp, yd, dy, ddy, &
                nr, itmax, toler, iter, intg, &
                nmat, materf, r, rini, epstr, &
                typess, essai, icomp, nvi, &
                vinf, &
                iret)
!     IRET = 0 CONVERGENCE
!          = 1 ITERATION SUIVANTE
!          = 2 RE-INTEGRATION
!          = 3 REDECOUPAGE DU PAS DE TEMPS
!
    if (iret .eq. 1) then
        goto 1
    else if (iret .eq. 2) then
        goto 2
    else if (iret .eq. 3) then
        goto 3
    end if
!
!     CONVERGENCE > INCREMENTATION DE  YF = YD + DY
    yf(1:ndt+nvi) = yd(1:ndt+nvi)+dy(1:ndt+nvi)
!
!     MISE A JOUR DE SIGF , VINF
    sigf(1:ndt) = yf(1:ndt)
!
!     POST-TRAITEMENTS POUR DES LOIS PARTICULIERES
    call lcplnf(BEHinteg, &
                rela_comp, vind, nbcomm, nmat, cpmono, &
                materf, iter, nvi, itmax, &
                toler, pgl, nfs, nsg, toutms, &
                hsr, dt, dy, yd, yf, &
                vinf, sigd, sigf, &
                deps, nr1, mod, timef, &
                iret)
!
    if (iret .ne. 0) then
        if (iret .eq. 2) goto 2
        goto 3
    end if
!
!     CONVERGENCE
    codret = 0
    goto 999
!
3   continue
!     NON CV, OU PB => REDECOUPAGE (LOCAL OU GLOBAL) DU PAS DE TEMPS
    codret = 1
!
999 continue
!
end subroutine
