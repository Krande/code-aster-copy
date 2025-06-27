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
!
#include "asterf_types.h"
!
interface
    subroutine asretm(lmasym, jtmp2, lgtmp2, nbterm, jsmhc,&
                      jsmdi, i1, i2)
        aster_logical :: lmasym
        integer(kind=8) :: jtmp2
        integer(kind=8) :: lgtmp2
        integer(kind=8) :: nbterm
        integer(kind=8) :: jsmhc
        integer(kind=8) :: jsmdi
        integer(kind=8) :: i1
        integer(kind=8) :: i2
    end subroutine asretm
end interface
