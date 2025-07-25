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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nonlinDSEnergyCreate(ds_energy)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/ismaem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/nonlinDSColumnVoid.h"
!
    type(NL_DS_Energy), intent(out) :: ds_energy
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Energy management
!
! Create energy management datastructure
!
! --------------------------------------------------------------------------------------------------
!
! Out ds_energy        : datastructure for energy management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_col
    integer(kind=8), parameter :: nb_col_defi = 8
    type(NL_DS_Table) :: table
    type(NL_DS_Column) :: column
!
    character(len=16), parameter :: cols_name(nb_col_defi) = (/ &
                                    'NUME_REUSE', 'INST      ', 'TRAV_EXT  ', &
                                    'ENER_CIN  ', 'ENER_TOT  ', 'TRAV_AMOR ', &
                                    'TRAV_LIAI ', 'DISS_SCH  '/)
!
    character(len=8), parameter :: cols_type(nb_col_defi) = (/ &
                                   'I  ', 'R  ', 'R  ', &
                                   'R  ', 'R  ', 'R  ', &
                                   'R  ', 'R  '/)
!
! --------------------------------------------------------------------------------------------------
!
    table%table_io%tablSymbName = 'PARA_CALC'
!
! - Create table
!
    table%nb_cols = nb_col_defi
    ASSERT(table%nb_cols .le. table%nb_cols_maxi)
    do i_col = 1, nb_col_defi
        call nonlinDSColumnVoid(column)
        column%name = cols_name(i_col)
        column%l_vale_inte = ASTER_FALSE
        column%l_vale_real = ASTER_FALSE
        column%vale_inte = ismaem()
        column%vale_real = r8vide()
        if (cols_type(i_col) .eq. 'I') then
            column%l_vale_inte = ASTER_TRUE
            column%vale_inte = 0
        elseif (cols_type(i_col) .eq. 'R') then
            column%l_vale_real = ASTER_TRUE
            column%vale_real = 0.d0
        else
            ASSERT(.false.)
        end if
        table%cols(i_col) = column
        table%l_cols_acti(i_col) = ASTER_FALSE
        table%indx_vale(i_col) = i_col
    end do
!
! - Set main parameters
!
    ds_energy%l_comp = ASTER_FALSE
    ds_energy%command = ' '
    ds_energy%table = table
!
end subroutine
