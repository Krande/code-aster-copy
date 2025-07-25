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

subroutine ComputeTableWidth(table, line_width, nb_cols_active)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Table), intent(in) :: table
    integer(kind=8), intent(out) :: line_width
    integer(kind=8), intent(out) :: nb_cols_active
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Table management
!
! Compute line of table
!
! --------------------------------------------------------------------------------------------------
!
! In  table            : datastructure for table
! Out line_width       : width of line
! Out nb_cols_active   : number of active columns
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_col, nb_cols
!
! --------------------------------------------------------------------------------------------------
!
    nb_cols = table%nb_cols
    line_width = 0
    nb_cols_active = 0
!
! - Number of active columns
!
    do i_col = 1, nb_cols
        if (table%l_cols_acti(i_col)) then
            nb_cols_active = nb_cols_active+1
        end if
    end do
!
! - Compute width
!
    line_width = 1
    do i_col = 1, nb_cols
        if (table%l_cols_acti(i_col)) then
            line_width = line_width+(16+1)
        end if
    end do
    ASSERT(line_width .le. 512)
!
end subroutine
