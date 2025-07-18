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
!
subroutine lcrolo(fami, kpg, ksp, mate, option, &
                  carcri, fm, df, vim, vip, &
                  taup, dtaudf, iret)
!
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterfort/gdclci.h"
#include "asterfort/gdclco.h"
#include "asterfort/gdclel.h"
#include "asterfort/gdclin.h"
#include "asterfort/gdcltg.h"
#include "asterfort/lcrofg.h"
#include "asterfort/lcrohy.h"
#include "asterfort/lcroma.h"
#include "asterfort/lcrotg.h"
#include "asterfort/lcroy1.h"
#include "asterfort/lcroy2.h"
#include "asterfort/lcroyi.h"
#include "blas/dcopy.h"
    integer(kind=8) :: kpg, ksp, mate, iret
    character(len=*) :: fami
    character(len=16) :: option
    real(kind=8) :: carcri(3)
    real(kind=8) :: fm(3, 3), df(3, 3), vim(9)
    real(kind=8) :: vip(9), taup(6), dtaudf(6, 3, 3)
!
!......................................................................
!       INTEGRATION DE LA LOI DE ROUSSELIER LOCAL
! EN GRANDES DEFORMATIONS DE TYPE NOUVELLE FORMULATION DE CANO-LORENTZ
!......................................................................
! IN  MATE    : ADRESSE DU MATERIAU CODE
! IN  OPTION  : OPTION DE CALCUL
! IN  CARCRI  : PARAM7TRES POUR L INTEGRATION DE LA LOI DE COMMPORTEMENT
!                CARCRI(1) = NOMBRE D ITERATIONS
!                CARCRI(3) = PRECISION SUR LA CONVERGENCE
! IN  FM      : GRADIENT DE LA TRANSFORMATION A L INSTANT PRECEDENT
! IN  DF      : INCREMENT DU GRADIENT DE LA TRANSFORMATION
! IN  VIM     : VARIABLES INTERNES A L INSTANT DU CALCUL PRECEDENT
!         VIM(1)   = P (DEFORMATION PLASTIQUE CUMULEE)
!         VIM(2)   = POROSITE
!         VIM(3)   = INDICATEUR DE PLASTICITE
!                  = 0 SOLUTION ELASTIQUE
!                  = 1 SOLUTION PLASTIQUE
!                  = 2 SOLUTION PLASTIQUE SINGULIERE
!         VIM(4:9) = DEFORMATION ELASTIQUE EULERIENNE EE = (ID-BE)/2)
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT TAUP    : CONTRAINTE A L INSTANT ACTUEL
! OUT DTAUDF  : DERIVEE DE TAU PAR RAPPORT A DF  * TRANSPOSE(DF)
! OUT IRET    : CODE RETOUR SUR L INTEGRATION DE LA LDC
!               SI LA FONCTION FONC = D*POROM*EXP(-K*TRETR/SIG1)
!               EST TROP GRANDE OU TROP PETITE ON REDECOUPE GLOBALEMENT
!               LE PAS DE TEMPS
! ----------------------------------------------------------------------
!  COMMON LOI DE COMPORTEMENT ROUSSELIER
!
    integer(kind=8) :: itemax, jprolp, jvalep, nbvalp
    real(kind=8) :: prec, young, nu, sigy, sig1, rousd, f0, fcr, acce
    real(kind=8) :: pm, rpm, fonc, fcd, dfcddj, dpmaxi, typoro
    common/lcrou/prec, young, nu, sigy, sig1, rousd, f0, fcr, acce,&
     &               pm, rpm, fonc, fcd, dfcddj, dpmaxi, typoro,&
     &               itemax, jprolp, jvalep, nbvalp
! ----------------------------------------------------------------------
!  COMMON GRANDES DEFORMATIONS CANO-LORENTZ
!
    integer(kind=8) :: ind1(6), ind2(6)
    real(kind=8) :: kr(6), rac2, rc(6)
    real(kind=8) :: lambda, mu, deuxmu, unk, troisk, cother
    real(kind=8) :: jm, dj, jp, djdf(3, 3)
    real(kind=8) :: etr(6), dvetr(6), eqetr, tretr, detrdf(6, 3, 3)
    real(kind=8) :: dtaude(6, 6)
!
    common/gdclc/&
     &          ind1, ind2, kr, rac2, rc,&
     &          lambda, mu, deuxmu, unk, troisk, cother,&
     &          jm, dj, jp, djdf,&
     &          etr, dvetr, eqetr, tretr, detrdf,&
     &          dtaude
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    character(len=9) :: typorot
    character(len=1) :: poum
    aster_logical :: resi, rigi, elas
    integer(kind=8) :: ij, indice
    real(kind=8) :: infini, petit
    real(kind=8) :: porom, poro, em(6), ep(6)
    real(kind=8) :: y, ym, x, seuil, dseuil, s, dp
    blas_int :: b_incx, b_incy, b_n
!
!    parameter (typoro = 'IMPLICITE')
!    parameter (typoro = 'EXPLICITE')
! ----------------------------------------------------------------------
!
!
! 1 - INITIALISATION
!
!    DONNEES DE CONTROLE DE L'ALGORITHME
    infini = r8gaem()
    petit = r8prem()
    itemax = nint(carcri(1))
    prec = carcri(3)
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    elas = option(11:14) .eq. 'ELAS'
    call gdclin()
!
!    LECTURE DES VARIABLES INTERNES
    pm = vim(1)
    porom = vim(2)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vim(4), b_incx, em, b_incy)
!
!    INITIALISATION SPECIFIQUE A RIGI_MECA_*
    if (resi) then
        poum = '+'
    else
        poum = '-'
        dp = 0.d0
        indice = nint(vim(3))
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, em, b_incx, ep, b_incy)
    end if
!
!    CARACTERISTIQUES MATERIAU
    call lcroma(fami, kpg, ksp, poum, mate)
!
!
    if (typoro .eq. 1) typorot = 'IMPLICITE'
    if (typoro .eq. 2) typorot = 'EXPLICITE'
!
    call gdclel(fami, kpg, ksp, poum, mate, &
                young, nu)
!
!
! 2 - CALCUL DES ELEMENTS CINEMATIQUES
!
    call gdclci(fm, df, em)
!
!
! 3 - CALCUL DE LA POROSITE ET TESTS ASSOCIES
!
    if (typorot .eq. 'EXPLICITE') then
!
!
!
!      POROSITE EXPLICITE
        porom = max(f0, porom)
        poro = porom
        fcd = rousd*poro
        dfcddj = 0.d0
!
    else if (typorot .eq. 'IMPLICITE') then
!
!    POROSITE FONCTION DE J
        poro = 1-(1-f0)/jp
        if (poro .gt. f0) then
            dfcddj = rousd*(1-f0)/jp**2
        else
            poro = f0
            dfcddj = 0
        end if
        fcd = rousd*poro
!
    end if
!
    if ((unk*abs(tretr)/sig1) .ge. 500.d0) then
        iret = 1
        goto 999
    end if
    fonc = fcd*exp(-unk*tretr/sig1)*exp(-cother/sig1)
!
    if ((fonc .ge. infini) .or. (fonc .le. petit)) then
        iret = 1
        goto 999
    end if
!
! 4 - INTEGRATION DE LA LOI DE COMPORTEMENT
!     PAR METHODE DE NEWTON AVEC BORNES CONTROLEES ET DICHOTOMIE
!  RESOLUTION DES EQATIONS:
!  - SI SEUIL(0)<0 => LA SOLUTION EST ELASTIQUE SINON
!  - SI S(0)>0     => LA SOLUTION EST PLASTIQUE ET REGULIERE
!                     ON RESOUD SEUIL(Y)=0
!  - SI S(0)<0     => ON RESOUD S(YS)=0
!                     YS EST SOLUTION SINGULIERE SI DP>2*EQ(DVE-DVETR)/3
!  - SINON ON RESOUD SEUIL(Y)=0 POUR Y>YS
!  AVEC  SEUIL(Y)= 2*MU*EQETR-S(Y)-3*MU*DP(Y)
!        Y       = K*X/SIG1
!        X       = TRE-TRETR
!        DP      = (Y*SIG1/K)*EXP(Y)/FONC
!        S(Y)    = -SIG1*FONC*EXP(-Y)+R(PM+DP)
!
    if (resi) then
!
! 4.1 - EXAMEN DE LA SOLUTION ELASTIQUE (Y=0)
!       LCROFG = CALCUL DU SEUIL ET DE SA DERIVEE
!                IN : Y - OUT : DP,S,SEUIL,DSEUIL
!
        y = 0
        call lcrofg(y, dp, s, seuil, dseuil)
        if (seuil .le. 0.d0) then
            indice = 0
            goto 600
        end if
!
! 4.2 - RESOLUTION SEUIL(Y)=0 QUAND S(0)>0
!       CALCUL DE Y PUIS DU DP CORRESPONDANT
!
        if (s .gt. 0) then
            y = lcroy1()
            call lcrofg(y, dp, s, seuil, dseuil)
            indice = 1
            goto 600
        end if
!
! 4.3 - EXAMEN DE LA SOLUTION SINGULIERE ( S(0)<0 )
! 4.3.1 - RESOLUTION S(Y)=0
!         CALCUL DE Y PUIS DU DP CORRESPONDANT
!
        y = lcroyi()
        call lcrofg(y, dp, s, seuil, dseuil)
!
! 4.3.2 - CONDITION POUR SOLUTION SINGULIERE
!
        if (2*eqetr/3.d0-dp .le. 0) then
            indice = 2
            goto 600
        end if
        ym = y
!
! 4.4 - RESOLUTION SEUIL(Y)=0 QUAND S(0)<0 : ON A S(YM)=0
!       CALCUL DE Y PUIS DU DP CORRESPONDANT
!
        y = lcroy2(ym)
        call lcrofg(y, dp, s, seuil, dseuil)
        indice = 1
!
600     continue
        x = sig1*y/unk
!
! 4.5 - CALCUL DE LA DEFORMATION ELASTIQUE
!
        if (indice .eq. 0) then
            b_n = to_blas_int(6)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, etr, b_incx, ep, b_incy)
        else
            do ij = 1, 6
                ep(ij) = (x+tretr)/3.d0*kr(ij)
            end do
            if (indice .eq. 1 .and. eqetr .gt. petit) then
                do ij = 1, 6
                    ep(ij) = ep(ij)+dvetr(ij)*(1.d0-3.d0*dp/(2.d0*eqetr))
                end do
            end if
        end if
!
! 4.6 - VALIDITE DE L'INTEGRATION
!
        if (dp .gt. dpmaxi) then
            iret = 2
        end if
!
! FIN DE LA RESOLUTION
!
!
! 5 - INTEGRATION DE LA POROSITE EN EXPLICITE
!
        if (typorot .eq. 'EXPLICITE') then
            poro = 1.d0-(1.d0-porom)*exp(-x)
            if (poro .ge. fcr) then
                poro = 1.d0-(1.d0-porom)*exp(-acce*x)
            end if
        end if
!
! 6 - CALCUL DES CONTRAINTES
!
        call gdclco(ep, taup)
!
! 7 - STOCKAGE DES VARIABLES INTERNES EN T+
!
        vip(1) = pm+dp
        vip(2) = poro
        vip(3) = indice
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, ep, b_incx, vip(4), b_incy)
!
! 8 - CORRECTION DE LA LOI D'ECOULEMENT A POSTERIORI
!
        if (indice .ge. 1) then
            call lcrohy(x, dp, em, vip(4))
        end if
!
    end if
!
! 9 - CALCUL DE LA MATRICE TANGENTE
!
    if (rigi) then
!
!      SI ON FORCE A UTILISER LA MATRICE DE DECHARGE
        if (elas) indice = 0
!
!      GRANDEURS DERIVEES COMMUNES A CANO-LORENTZ
        call gdcltg(df, ep)
!
!      DERIVATIONS SPECIFIQUES A LA LOI DE ROUSSELIER
        call lcrotg(indice, dp, ep, dtaudf)
!
    end if
!
999 continue
end subroutine
