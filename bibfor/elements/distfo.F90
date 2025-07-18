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
function distfo(zimat, kfonc, xx, yy, normx, &
                normy)
!
    implicit none
!
!     CALCUL
!
! IN  ZIMAT :
! IN  KFONC :
! IN  XX :
! IN  YY :
!
! OUT NORMX :
! OUT NORMY :
!
#include "asterfort/cdnfon.h"
#include "asterfort/rcvalb.h"
    integer(kind=8) :: i, itmax, ier, zimat, kpg, spt
!
    real(kind=8) :: distfo, xx, yy, normx, normy, x0, y0
    real(kind=8) :: xi, yi, xm1, res, ym1, rp, dym1, tol
    real(kind=8) :: xm2, ym2, dyi, val(1)
!
    integer(kind=8) :: codres(1)
    character(len=8) :: kfonc, fami, poum
    character(len=16) :: phenom
!
    phenom = 'GLRC_DAMAGE'
!
    tol = 1.0d-3
    x0 = xx/normx
    y0 = yy/normy
    res = 1.0d6
    xi = x0
!
    ym1 = 0.0d0
    xi = 0.0d0
    itmax = 1000
    xm1 = 1.0d20
    ym1 = 1.0d20
!
    xi = xi*normx
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zimat, &
                ' ', phenom, 1, 'X ', [xi], &
                1, kfonc, val, codres, 1)
    yi = val(1)
    call cdnfon(zimat, kfonc, xi, 1, dyi, &
                ier)
    yi = yi/normy
    xi = xi/normx
    dyi = dyi*normx/normy
!
    do i = 1, itmax
        xm2 = xm1
        ym2 = ym1
        xm1 = xi
        ym1 = yi
        dym1 = dyi
!
        rp = (xm2-xm1)*(xm2-xm1)+(ym2-ym1)*(ym2-ym1)
        res = sqrt(rp*rp)
!
        if (res .lt. tol) goto 30
!
        rp = dym1/(dym1*dym1+1.0d0)
        xi = rp*(y0-ym1+dym1*xm1)+x0/(dym1*dym1+1.0d0)
        xi = xi*normx
!
        call rcvalb(fami, kpg, spt, poum, zimat, &
                    ' ', phenom, 1, 'X ', [xi], &
                    1, kfonc, val, codres, 1)
        yi = val(1)
        call cdnfon(zimat, kfonc, xi, 1, dyi, &
                    ier)
!
        yi = yi/normy
        xi = xi/normx
        dyi = dyi*normx/normy
!
    end do
!
30  continue
!
    rp = (xm1-x0)*(xm1-x0)
    rp = rp+(ym1-y0)*(ym1-y0)
    distfo = sqrt(rp)
!
end function
