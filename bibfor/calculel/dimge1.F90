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
function dimge1(ige1, ige2)
    implicit none
    integer(kind=8) :: dimge1
#include "asterfort/assert.h"
    integer(kind=8) :: ige1, ige2
!     -- SERT AU CALCUL DE DIM_GEOM D'UN LIGREL
!    IN:
!       IGE1 : VALEUR DE DIM_GEOM (/1/2/3/120/023/103/123)
!       IGE2 : VALEUR DE DIM_GEOM (/1/2/3/120/023/103/123)
!    OUT:
!       DIMGE1   : "CUMUL" DE IGE1 ET IGE2
!
! ----------------------------------------------------------------------
    integer(kind=8) :: i1(3), i2(3), i3(3), k
! DEB ------------------------------------------------------------------
!
!     -- SI IGE1 (OU IGE2) = 0, C'EST FACILE :
    if (ige1 .eq. 0) then
        dimge1 = ige2
        goto 999
    end if
    if (ige2 .eq. 0) then
        dimge1 = ige1
        goto 999
    end if
!
!     -- ON DECODE IGE1 DANS I1 :
    do k = 1, 3
        i1(k) = 0
    end do
    if (ige1 .gt. 10) then
        i1(1) = ige1/100
        i1(3) = mod(ige1, 10)
        i1(2) = (mod(ige1, 100)-i1(3))/10
    else
        ASSERT(ige1 .ge. 0 .and. ige1 .le. 3)
        if (ige1 .eq. 1) i1(1) = 1
        if (ige1 .eq. 2) i1(2) = 2
        if (ige1 .eq. 3) i1(3) = 3
    end if
!
!     -- ON DECODE IGE2 DANS I2 :
    do k = 1, 3
        i2(k) = 0
    end do
    if (ige2 .gt. 10) then
        i2(1) = ige2/100
        i2(3) = mod(ige2, 10)
        i2(2) = (mod(ige2, 100)-i2(3))/10
    else
        ASSERT(ige2 .ge. 0 .and. ige2 .le. 3)
        if (ige2 .eq. 1) i2(1) = 1
        if (ige2 .eq. 2) i2(2) = 2
        if (ige2 .eq. 3) i2(3) = 3
    end if
!
!     -- ON CALCULE I3 :
    do k = 1, 3
        ASSERT(i1(k) .eq. 0 .or. i1(k) .eq. k)
        ASSERT(i2(k) .eq. 0 .or. i2(k) .eq. k)
        i3(k) = max(i1(k), i2(k))
    end do
!
!     -- ON RECODE LE RESULTAT :
    if ((i3(1)+i3(2)+i3(3)) .gt. 1) then
        dimge1 = 100*i3(1)
        dimge1 = dimge1+10*i3(2)
        dimge1 = dimge1+1*i3(3)
    else
        if (i3(1) .eq. 1) dimge1 = 1
        if (i3(2) .eq. 1) dimge1 = 2
        if (i3(3) .eq. 1) dimge1 = 3
    end if
!
999 continue
!
end function
