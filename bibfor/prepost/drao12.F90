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

subroutine drao12(coord1, coord2, xo1o2, yo1o2, zo1o2, &
                  do1o2, coorda, ray)
!
    implicit none
!
! DECLARATION GLOBALE
!
#include "asterfort/assert.h"
    real(kind=8) :: coord1(3)
    real(kind=8) :: coord2(3)
    real(kind=8) :: xo1o2, yo1o2, zo1o2, do1o2
    real(kind=8) :: coorda(3)
    real(kind=8) :: ray(2)
!
! DECLARATION LOCALE
!
    real(kind=8) :: xao1, yao1, zao1, xao2, yao2, zao2
    real(kind=8) :: dao1, dao2, pao1o2, pao2o1, n2, n3, n4
    real(kind=8) :: r, cos1, cos2
!
    xao1 = coorda(1)-coord1(1)
    yao1 = coorda(2)-coord1(2)
    zao1 = coorda(3)-coord1(3)
    xao2 = coorda(1)-coord2(1)
    yao2 = coorda(2)-coord2(2)
    zao2 = coorda(3)-coord2(3)
!
    dao1 = sqrt(xao1**2+yao1**2+zao1**2)
    dao2 = sqrt(xao2**2+yao2**2+zao2**2)
!
    r = sqrt( &
        (yao1*zo1o2-zao1*yo1o2)**2+(zao1*xo1o2-xao1*zo1o2)**2+(xao1*yo1o2-yao1*xo1&
        &o2)**2 &
        )/do1o2
!
    ray(1) = r
    ray(2) = 0.0d0
!
    pao1o2 = (xao1*xo1o2)+(yao1*yo1o2)+(zao1*zo1o2)
    if (dao1 .eq. 0.d0) then
        cos1 = 1.0d0
    else
        cos1 = pao1o2/(do1o2*dao1)
    end if
    pao2o1 = (-xao2*xo1o2)+(-yao2*yo1o2)+(-zao2*zo1o2)
    if (dao2 .eq. 0.d0) then
        cos2 = 1.0d0
    else
        cos2 = pao2o1/(do1o2*dao2)
    end if
    n2 = (dao1*cos1)/do1o2
    n3 = (dao2*cos2)/do1o2
    n4 = n2+n3
    if (abs(1.d0-n2) .le. 1.0d-10) n2 = 1.0d0
    if (abs(1.d0-n3) .le. 1.0d-10) n3 = 1.0d0
    if (abs(1.d0-n4) .le. 1.0d-10) n4 = 1.0d0
    ASSERT(n4 .eq. 1.0d0)
!
    if (n2 .gt. 1.0d0 .or. n3 .gt. 1.0d0) then
        if (n2 .gt. 1.0d0) ray(2) = -1.0d0
        if (n3 .gt. 1.0d0) ray(2) = 1.0d0
    end if
!
end subroutine
