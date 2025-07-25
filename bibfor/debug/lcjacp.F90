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

subroutine lcjacp(fami, kpg, ksp, rela_comp, toler, &
                  itmax, mod, imat, nmat, materd, &
                  materf, nr, nvi, timed, timef, &
                  deps, epsd, vind, vinf, yd, &
                  nbcomm, cpmono, pgl, nfs, &
                  nsg, toutms, hsr, dy, r, &
                  drdy, verjac, drdyb, iret, crit)
!
!     CONSTRUCTION DE LA MATRICE JACOBIENNE PAR PERTURBATION
!     IN  FAMI   :  FAMILLE DE POINT DE GAUSS
!         KPG    :  NUMERO DU POINT DE GAUSS
!         KSP    :  NUMERO DU SOUS-POINT DE GAUSS
!         LOI    :  MODELE DE COMPORTEMENT
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
!     VAR DEPS   :  INCREMENT DE DEFORMATION
!     IN  EPSD   :  DEFORMATION A T
!         SIGD   :  CONTRAINTE A T
!         VIND   :  VARIABLES INTERNES A T
!         VINF   :  VARIABLES INTERNES A T+DT
!         YD     :  VARIABLES A T   = ( SIGD  VIND  (EPSD3)   )
!         COMP   :  COMPORTEMENT
!         DY     :  INCREMENT DES VARIABLES = ( DSIG  DVIN  (DEPS3)  )
!         R      :  VECTEUR RESIDU
!         DRDY   :  JACOBIEN
!
!         VERJAC : =0 : PAS DE VERIFICATION
!         =1 : CONSTRUCTION DE LA JACOBIENNE PAR PERTURBATION (LCJACP)
!                COMPARAISON A LA MATRICE JACOBIENNE ISSU DE LCJACB
!         =2 : UTILISATION DE LA JACOBIENNE PAR PERTURBATION (LCJACP)
!                COMME MATRICE JACOBIENNE A LA PLACE DE LCJACB
!     OUT DRDYB  : MATRICE JACOBIENNE PAR PERTURBATION
! ----------------------------------------------------------------------
! aslint: disable=W1306,W1504
    implicit none
!
#include "asterc/r8miem.h"
#include "asterfort/lcresi.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nmat, nbcomm(nmat, 3), nr, impr, vali(2), nfs, nsg
    integer(kind=8) :: imat, i, j, itmax, iret, kpg, ksp, nvi, verjac
!
    real(kind=8) :: toler, epsd(6), deps(6), vind(nvi), vinf(nvi), timed, timef
    real(kind=8) :: err
!
!     DIMENSIONNEMENT DYNAMIQUE (MERCI F90)
    real(kind=8) :: dy(nr), r(nr), drdyb(nr, nr), rini(nr), dyini(nr), rp(nr)
    real(kind=8) :: rm(nr)
    real(kind=8) :: drdy(nr, nr), yd(nr), dym(nr), dyp(nr), yfp(nr), yfm(nr)
!
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), pgl(3, 3), eps1, eps2
    real(kind=8) :: eps0
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg), crit(*)
    real(kind=8) :: valr(4), maxtgt, normd1, normd2, maxerr
!
    character(len=8) :: mod
    character(len=16) :: rela_comp
    character(len=24) :: cpmono(5*nmat+1)
    character(len=*) :: fami
    data impr/0/
! ----------------------------------------------------------------------
    dyini(1:nr) = dy(1:nr)
    rini(1:nr) = r(1:nr)
    maxtgt = 0.d0
    normd1 = 0.d0
    normd2 = 0.d0
!
    do i = 1, 6
        normd1 = normd1+dyini(i)*dyini(i)
    end do
!
    do i = 7, nr
        normd2 = normd2+dyini(i)*dyini(i)
    end do
!
    if (normd1 .lt. r8miem()) then
        do i = 1, 6
            normd1 = normd1+yd(i)*yd(i)
        end do
    end if
    if (normd2 .lt. r8miem()) then
        do i = 7, nr
            normd2 = normd2+yd(i)*yd(i)
        end do
    end if
!
    eps0 = 1.d-7
    eps1 = eps0
    eps2 = eps0
    if (normd1 .gt. r8miem()) then
        eps1 = eps1*sqrt(normd1)
    end if
    if (normd2 .gt. r8miem()) then
        eps2 = eps2*sqrt(normd2)
    end if
!
    do i = 1, nr
        dyp(1:nr) = dyini(1:nr)
        if (i .le. 6) then
            dyp(i) = dyp(i)+eps1
        else
            dyp(i) = dyp(i)+eps2
        end if
        yfp(1:nr) = yd(1:nr)+dyp(1:nr)
        call lcresi(fami, kpg, ksp, rela_comp, mod, &
                    imat, nmat, materd, materf, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, nr, nvi, vind, &
                    vinf, itmax, toler, timed, timef, &
                    yd, yfp, deps, epsd, dyp, &
                    rp, iret, crit)
        if (iret .gt. 0) then
            goto 999
        end if
        dym(1:nr) = dyini(1:nr)
        if (i .le. 6) then
            dym(i) = dym(i)-eps1
        else
            dym(i) = dym(i)-eps2
        end if
        yfm(1:nr) = yd(1:nr)+dym(1:nr)
        call lcresi(fami, kpg, ksp, rela_comp, mod, &
                    imat, nmat, materd, materf, &
                    nbcomm, cpmono, pgl, nfs, nsg, &
                    toutms, hsr, nr, nvi, vind, &
                    vinf, itmax, toler, timed, timef, &
                    yd, yfm, deps, epsd, dym, &
                    rm, iret, crit)
        if (iret .gt. 0) then
            goto 999
        end if
!        SIGNE - CAR LCRESI CALCULE -R
        do j = 1, nr
            if (i .le. 6) then
                drdyb(j, i) = -(rp(j)-rm(j))/2.d0/eps1
            else
                drdyb(j, i) = -(rp(j)-rm(j))/2.d0/eps2
            end if
        end do
    end do
!
! COMPARAISON DRDY ET DRDYB
!
    maxerr = 0.d0
    err = 0.d0
    if ((verjac .eq. 1) .and. (impr .eq. 0)) then
        do i = 1, nr
            do j = 1, nr
                if (abs(drdy(i, j)) .gt. maxtgt) then
                    maxtgt = abs(drdy(i, j))
                end if
            end do
        end do
        do i = 1, nr
            do j = 1, nr
                if (abs(drdy(i, j)) .gt. (1.d-9*maxtgt)) then
                    if (abs(drdyb(i, j)) .gt. (1.d-9*maxtgt)) then
                        err = abs(drdy(i, j)-drdyb(i, j))/drdyb(i, j)
                        if (err .gt. 1.d-3) then
                            vali(1) = i
                            vali(2) = j
!
                            valr(1) = timef
                            valr(2) = err
                            valr(3) = drdyb(i, j)
                            valr(4) = drdy(i, j)
                            call utmess('I', 'DEBUG_1', ni=2, vali=vali, nr=4, &
                                        valr=valr)
                            maxerr = max(maxerr, abs(err))
                            impr = 1
                        end if
                    end if
                end if
            end do
        end do
    end if
!
!     UTILISATION DE DRDYB COMME MATRICE JACOBIENNE
    if (verjac .eq. 2) then
        drdy(:, :) = drdyb(:, :)
    end if
!
999 continue
end subroutine
