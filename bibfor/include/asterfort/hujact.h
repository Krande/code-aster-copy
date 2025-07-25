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
    subroutine hujact(mater, vind, vinf, vins, sigd,&
                      sigf, negmul, chgmec, indi)
        real(kind=8) :: mater(22, 2)
        real(kind=8) :: vind(*)
        real(kind=8) :: vinf(*)
        real(kind=8) :: vins(50)
        real(kind=8) :: sigd(6)
        real(kind=8) :: sigf(6)
        aster_logical :: negmul(8)
        aster_logical :: chgmec
        integer(kind=8) :: indi(7)
    end subroutine hujact
end interface
