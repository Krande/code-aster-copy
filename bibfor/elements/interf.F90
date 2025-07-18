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
subroutine interf(mater, kfonc1, kfonc2, normf, x0, &
                  xrac)
    implicit none
#include "asterfort/cdnfo2.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
    character(len=16) :: kfonc1, kfonc2
    character(len=8) :: mater, k8b
    real(kind=8) :: normf, x0
    real(kind=8) :: xrac
    integer(kind=8) :: ier1, ier2, iter, itermx
    real(kind=8) :: fx1(1), fx2(1), fx, dfx1, dfx2, dfx, tole, err
    character(len=32) :: phenom
    integer(kind=8) :: icodr2(1)
!
    phenom = 'GLRC_DAMAGE'
    iter = 0
    itermx = 100
    xrac = x0
!
    k8b = 'X '
    call rcvale(mater, phenom, 1, k8b, [xrac], &
                1, kfonc1, fx1(1), icodr2(1), 1)
    call rcvale(mater, phenom, 1, k8b, [xrac], &
                1, kfonc2, fx2(1), icodr2(1), 1)
!
    fx = fx1(1)-fx2(1)
    err = abs(fx)
    tole = 1.d-8*normf
!
    do iter = 1, itermx
        if (err .le. tole) goto 10
        call cdnfo2(mater, kfonc1, xrac, 1, dfx1, &
                    ier1)
        call cdnfo2(mater, kfonc2, xrac, 1, dfx2, &
                    ier2)
        dfx = dfx1-dfx2
!
        if ((abs(dfx) .lt. 1.d-12) .or. (ier1 .gt. 0) .or. (ier2 .gt. 0)) then
            call utmess('F', 'ELEMENTS2_27')
        end if
!
        xrac = xrac-fx/dfx
        call rcvale(mater, phenom, 1, k8b, [xrac], &
                    1, kfonc1, fx1(1), icodr2(1), 1)
        call rcvale(mater, phenom, 1, k8b, [xrac], &
                    1, kfonc2, fx2(1), icodr2(1), 1)
        fx = fx1(1)-fx2(1)
        err = abs(fx)
    end do
    call utmess('F', 'ELEMENTS2_27')
!
10  continue
!
end subroutine
