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
#include "asterf_types.h"
!
interface
    subroutine ornorm(mesh, listCellNume, nbCell, reorie, norien, nconex, &
                      onlySkin1D_)
        character(len=8), intent(in) :: mesh
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), pointer :: listCellNume(:)
        aster_logical, intent(in) :: reorie
        integer(kind=8), intent(out) :: norien, nconex
        aster_logical, intent(in), optional :: onlySkin1D_
    end subroutine ornorm
end interface
