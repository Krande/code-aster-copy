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
#include "MeshTypes_type.h"
!
interface
    subroutine pjxxut(projDime, typeSelect, &
                      entity1, entity2, &
                      nbCellSelect1, listCellSelect1, &
                      nbNodeSelect2, listNodeSelect2, &
                      mesh1, mesh2, &
                      nbCellType, cellListNume, cellListCode)
        character(len=2), intent(in) :: projDime
        character(len=*), intent(in) :: typeSelect
        character(len=8), intent(in) :: entity1, entity2
        integer(kind=8), intent(in) :: nbCellSelect1, listCellSelect1(*)
        integer(kind=8), intent(in) :: nbNodeSelect2, listNodeSelect2(*)
        character(len=8), intent(out) :: mesh1, mesh2
        integer(kind=8), intent(out) :: nbCellType, cellListNume(MT_NTYMAX)
        character(len=8), intent(out) :: cellListCode(MT_NTYMAX)
    end subroutine pjxxut
end interface
