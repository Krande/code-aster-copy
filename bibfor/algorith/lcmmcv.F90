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
subroutine lcmmcv(yd, dy, ddy, nr, itmax, &
                  toler, iter, r, rini, epstr, &
                  irteti)
! person_in_charge: jean-michel.proix at edf.fr
! aslint: disable=W1306
    implicit none
!       MONOCRISTAL        : CONTROLE DE LA CONVERGENCE
!                                  DE LA CONFORMITE DE LA SOLUTION DP
!                                  ET DE LA RE-INTEGRATION
!                                  ET DU REDECOUPAGE DU PAS DE TEMPS
!                                  SUR LA NORME DU RESIDU
!       ----------------------------------------------------------------
!       IN
!            DY     :  VECTEUR SOLUTION DY = ( DSIG DVINT)
!            DDY    :  VECTEUR CORRECTION SUR LA SOLUTION
!            NR     :  DIMENSION DY DDY
!            ITMAX  :  NB MAXI D ITERATIONS LOCALES
!            TOLER  :  TOLERANCE A CONVERGENCE
!            ITER   :  NUMERO ITERATION COURANTE
!            R      :  R(Y) RESIDU A L'ITERATION COURANTE
!            RINI   :  R(Y0) RESIDU A L'ITERATION 1
!       OUT  IRTETI  :  =0 CONVERGENCE
!                       =1 ITERATIONS SUPPLEMENTAIRE (ITER<ITMAX)
!                       =3 ITMAX ATTEINT REDECOUPAGE
!       ----------------------------------------------------------------
#include "asterc/r8prem.h"
    integer(kind=8) :: itmax, iter, nr, irteti, ndt, ndi, i
    real(kind=8) :: toler, yd(nr), ddy(nr), dy(nr), r(nr), rini(nr), epstr(6)
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    real(kind=8) :: errdy(nr), errr(nr), e1, e2, e1ini, e2ini
!       ----------------------------------------------------------------
!
! -   EVALUATION  DE L'ERREUR  EN RESIDU (DEFORMATIONS)
!
!     TYP  : TYPE D'ERREUR A CALCULER
!          0 =  MAX | DDY /DY |     < EPS
!          1 = || DDY || / || DY || < EPS non si termes d'ordres /=
!          2 = || DDYi / DYi ||     < EPS
!      CALL LCVERR ( RINI, R, 6, 1, ERRR  )
!
    e1 = 0.d0
    e2 = 0.d0
    e1ini = 0.d0
    e2ini = 0.d0
    do i = 1, 6
        e1 = max(e1, abs(r(i)))
!         E1INI = MAX(E1INI, ABS(RINI(I)))
        e1ini = max(e1ini, abs(epstr(i)))
    end do
!
    errr(1) = e1
    if (e1ini .gt. r8prem()) then
        errr(1) = e1/e1ini
    end if
!
    do i = 7, nr
        e2 = max(e2, abs(r(i)))
!         E2INI = MAX(E2INI, ABS(RINI(I)))
        e2ini = max(e2ini, abs(yd(i)+dy(i)))
    end do
!
    errr(2) = e2
    if (e2ini .gt. r8prem()) then
        errr(2) = e2/e2ini
    end if
!
!     MAX DES 6 PREMIERS TERMES ET DES SUIVANTS
    errr(1) = max(errr(1), errr(2))
!
!     ERREUR SUR LA SOLUTION
!
    errdy(1) = 0.d0
    errdy(3) = 0.d0
    do i = 1, 6
        errdy(1) = max(errdy(1), abs(ddy(i)))
        errdy(3) = max(errdy(3), abs(dy(i)))
    end do
    if (errdy(3) .gt. r8prem()) then
        errdy(1) = errdy(1)/errdy(3)
    end if
    errdy(2) = 0.d0
    errdy(4) = 0.d0
    do i = ndt+1, nr
        errdy(2) = max(errdy(2), abs(ddy(i)))
        errdy(4) = max(errdy(3), abs(dy(i)))
    end do
    if (errdy(4) .gt. r8prem()) then
        errdy(2) = errdy(2)/errdy(4)
    end if
!
!     MAX DES 6 PREMIERS TERMES ET DES SUIVANTS
    errdy(1) = max(errdy(1), errdy(2))
!
!
!
    if (iter .le. itmax) then
! -             CONVERGENCE
        if (errr(1) .le. toler) then
            irteti = 0
            goto 999
        end if
        if (errdy(1) .le. toler) then
            irteti = 0
            goto 999
        end if
!
! -     NON CONVERGENCE ITERATION SUIVANTE
        irteti = 1
        goto 999
    else
! -     NB ITERATION MAXIMUM ATTEINT SANS CONVERGENCE
        irteti = 3
    end if
!
999 continue
end subroutine
