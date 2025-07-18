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
subroutine lcreli(fami, kpg, ksp, rela_comp, mod, &
                  imat, nmat, materd, materf, nbcomm, &
                  cpmono, pgl, nfs, nsg, toutms, &
                  hsr, nr, nvi, vind, vinf, &
                  itmax, toler, timed, timef, yd, &
                  yf, deps, epsd, dy, r, &
                  ddy, iret, crit)
! aslint: disable=W1306,W1504
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/lcresi.h"
#include "blas/ddot.h"
    real(kind=8) :: ddy(*)
!
!     VARIABLES EN ARGUMENT DE LCRESI
    integer(kind=8) :: imat, nmat, nr, nvi, kpg, ksp, itmax, nfs, nsg
    real(kind=8) :: deps(6), epsd(6), vind(*), toler, vinf(*)
    real(kind=8) :: r(*), yd(*), yf(*), dy(*)
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: timed, timef
    character(len=8) :: mod
    character(len=16) :: rela_comp
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg), crit(*)
    character(len=*) :: fami
!
    integer(kind=8) :: nbcomm(nmat, 3)
    real(kind=8) :: pgl(3, 3)
    character(len=24) :: cpmono(5*nmat+1)
!
!
!
!     ----------------------------------------------------------------
!     RECHERCHE LINEAIRE POUR PLASTI
!
!     ON CREE UNE FONCTIONNELLE F(X) = 1/2 || R(X) ||^2
!     ET ON CHERCHE RHO = ARGMIN (        F(DY+RHO.DDY)       )
!     CAD           RHO = ARGMIN ( 1/2 || R(DY+RHO.DDY) ||^2  )
!     ET ON MET A JOUR YF, R ET DY
!
!     ON UTILISE L'ALGORITHME AVEC REBROUSSEMENT ET LA REGLE D'ARMIJO
!     (W EST LA PARAMETRE DE LA REGLE D'ARMIJO)
!     AVEC RABATEMMENT SUR [RHOMIN,RHOMAX]
!     ON COMMENCE PAR UN ESSAI AVEC UN RHO = 1
!     PUIS UN ESSAI AVEC UNE INTERPOLATION QUADRATIQUE
!     PUIS PLUSIEURS ESSAIS (NB ESSAIS = IMXRHO) D'INTERPOLATION CUBIQUE
!
!     IN  TOUS LES ARGUMENTS DE LCRESI.F
!     IN  DDY    :  CORRECTION DE L'INCREMENT = DIRECTION DE DESCENTE
!     OUT R      :  VECTEUR RESIDU
!     OUT DY     :  INCREMENT DES VARIABLES
!     OUT YF     :  VARIABLES A T+DT
!     IN /OUT IRET : CODE RETOUR D'ERREUR (DIFFERENT DE 0 SI PB)
!
!     ----------------------------------------------------------------
    integer(kind=8) :: i, iret, itrho, imxrho
    real(kind=8) :: f, df, w, rhomin, rhomax
    real(kind=8) :: rhoddy(nr), dyp(nr), rp(nr), yfp(nr)
    real(kind=8) :: rho0, fp0, rho1, fp1, fp2, rho2, rho05, fsup
    real(kind=8) :: m(2, 2), s(2), a, b
    blas_int :: b_incx, b_incy, b_n
    parameter(w=1.d-4)
    parameter(rhomin=0.1d0, rhomax=0.5d0)
    parameter(imxrho=2)
!     ----------------------------------------------------------------
!
!     REMARQUE : ON POURRAIT METTRE DANS UNE ROUTINE UTILITAIRE LES 8
!     LIGNES CORRESPONDANTES AU CALCUL DU R ACTUALISE, MAIS CA VAUT PAS
!     VRAIMENT LE COUP
!
!
!     FONCTIONNELLE EN "MOINS" : F = 1/2 || R(DY) ||^2
    b_n = to_blas_int(nr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    f = 0.5d0*ddot(b_n, r, b_incx, r, b_incy)
!
!     DERIVEE DE LA FONCTIONNELLE EN "MOINS" : DF=<GRAD(F).DDY>=<R.DDY>
    b_n = to_blas_int(nr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    df = -ddot(b_n, r, b_incx, r, b_incy)
!
!     ------------------------------------
!     ESSAI AVEC LE PAS DE NEWTON RHO0 = 1
!     ------------------------------------
!
    rho0 = 1
!     CALCUL DE DY "PLUS" : DYP
    rhoddy = rho0*ddy(1:nr)
    dyp = rhoddy+dy(1:nr)
    yfp = yd(1:nr)+dyp
!     CALCUL DE R "PLUS" : RP
    call lcresi(fami, kpg, ksp, rela_comp, mod, &
                imat, nmat, materd, materf, nbcomm, &
                cpmono, pgl, nfs, nsg, toutms, &
                hsr, nr, nvi, vind, vinf, &
                itmax, toler, timed, timef, yd, &
                yfp, deps, epsd, dyp, rp, &
                iret, crit)
!
    if (iret .ne. 0) goto 999
!
!     TEST DE LA REGLE D'ARMIJO : SI TEST REUSSI, ON SORT
    b_n = to_blas_int(nr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    fp0 = 0.5d0*ddot(b_n, rp, b_incx, rp, b_incy)
    if (fp0 .lt. r8prem()) goto 888
    if (fp0 .le. f+w*rho0*df) goto 888
!
!     ------------------------------------
!     TEST SUPPLEMENTAIRE AVEC RHO = 0.5
!     ------------------------------------
!
    rho05 = 0.5d0
    rhoddy = rho05*ddy(1:nr)
    dyp = rhoddy+dy(1:nr)
    yfp = yd(1:nr)+dyp
    call lcresi(fami, kpg, ksp, rela_comp, mod, &
                imat, nmat, materd, materf, nbcomm, &
                cpmono, pgl, nfs, nsg, toutms, &
                hsr, nr, nvi, vind, vinf, &
                itmax, toler, timed, timef, yd, &
                yfp, deps, epsd, dyp, rp, &
                iret, crit)
    if (iret .ne. 0) goto 999
!
!     TEST DE LA REGLE D'ARMIJO : SI TEST REUSSI, ON SORT
    b_n = to_blas_int(nr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    fsup = 0.5d0*ddot(b_n, rp, b_incx, rp, b_incy)
    if (fsup .lt. r8prem()) goto 888
    if (fsup .le. f+w*0.5d0*df) goto 888
!
!     ----------------------------------------
!     INTERPOLATION QUADRATIQUE (ENTRE 0 ET 1)
!     ----------------------------------------
!
    ASSERT(abs(fp0-f-df) .gt. r8prem())
    rho1 = -0.5d0*df/(fp0-f-df)
!
!     PROJECTION SUR L'INTERVALLE [RHOMIN,RHOMAX]
    if (rho1 .lt. rhomin*rho0) rho1 = rhomin*rho0
    if (rho1 .gt. rhomax*rho0) rho1 = rhomax*rho0
!
    rhoddy = rho1*ddy(1:nr)
    dyp = rhoddy+dy(1:nr)
    yfp = yd(1:nr)+dyp
    call lcresi(fami, kpg, ksp, rela_comp, mod, &
                imat, nmat, materd, materf, nbcomm, &
                cpmono, pgl, nfs, nsg, toutms, &
                hsr, nr, nvi, vind, vinf, &
                itmax, toler, timed, timef, yd, &
                yfp, deps, epsd, dyp, rp, &
                iret, crit)
    if (iret .ne. 0) goto 999
!
!     TEST DE LA REGLE D'ARMIJO : SI TEST REUSSI, ON SORT
    b_n = to_blas_int(nr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    fp1 = 0.5d0*ddot(b_n, rp, b_incx, rp, b_incy)
    if (fp1 .lt. r8prem()) goto 888
    if (fp1 .le. f+w*rho1*df) goto 888
!
!     ------------------------------------
!     INTERPOLATIONS CUBIQUES
!     ------------------------------------
!
    do itrho = 1, imxrho
        m(1, 1) = 1.d0/(rho0**2)
        m(1, 2) = -1.d0/(rho1**2)
        m(2, 1) = -rho1/(rho0**2)
        m(2, 2) = rho0/(rho1**2)
        s(1) = fp0-f-df*rho0
        s(2) = fp1-f-df*rho1
        ASSERT(abs(rho0-rho1) .gt. r8prem())
        a = 1.d0/(rho0-rho1)*(m(1, 1)*s(1)+m(1, 2)*s(2))
        b = 1.d0/(rho0-rho1)*(m(2, 1)*s(2)+m(2, 2)*s(2))
        if (abs(3.d0*a) .le. r8prem()) goto 888
        rho2 = (-b+sqrt(b**2-3.d0*a*df))/(3.d0*a)
!
!       PROJECTION SUR L'INTERVALLE [RHOMIN,RHOMAX]
        if (rho2 .lt. rhomin*rho1) rho2 = rhomin*rho1
        if (rho2 .gt. rhomax*rho1) rho2 = rhomax*rho1
!
        rhoddy = rho2*ddy(1:nr)
        dyp = rhoddy+dy(1:nr)
        yfp = yd(1:nr)+dyp
        call lcresi(fami, kpg, ksp, rela_comp, mod, &
                    imat, nmat, materd, materf, nbcomm, &
                    cpmono, pgl, nfs, nsg, toutms, &
                    hsr, nr, nvi, vind, vinf, &
                    itmax, toler, timed, timef, yd, &
                    yfp, deps, epsd, dyp, rp, &
                    iret, crit)
        if (iret .ne. 0) goto 999
!
!       TEST DE LA REGLE D'ARMIJO : SI TEST REUSSI, ON SORT
        b_n = to_blas_int(nr)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        fp2 = 0.5d0*ddot(b_n, rp, b_incx, rp, b_incy)
        if (fp2 .lt. r8prem()) goto 888
        if (fp2 .le. f+w*rho2*df) goto 888
!
!       NOUVELLE INTERPOLATION CUBIQUE AVEC LES DEUX DERNIERS RHO
        rho0 = rho1
        rho1 = rho2
        fp0 = fp1
        fp1 = fp2
    end do
!
!     ON A FAIT TOUTES LES INTERATIONS D'INTERPOLATIONS CUBIQUES
!
888 continue
!
!     EN ECRASE LES ENTREES AVEC LES VARIABLES RE-ACTUALISEES
    do i = 1, nr
        r(i) = rp(i)
        yf(i) = yfp(i)
        dy(i) = dyp(i)
    end do
!
999 continue
end subroutine
