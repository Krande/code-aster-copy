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
subroutine nonlinDSEnergyInit(resultName, ds_energy)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/nonlinDSTableIOCreate.h"
#include "asterfort/nonlinDSTableIOSetPara.h"
#include "asterfort/nonlinDSTableIOGetName.h"
!
    character(len=8), intent(in) :: resultName
    type(NL_DS_Energy), intent(inout) :: ds_energy
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Energy management
!
! Initializations for energy
!
! --------------------------------------------------------------------------------------------------
!
! In  resultName       : name of results datastructure
! IO  ds_energy        : datastructure for energy management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_col
    character(len=24) :: col_name
    type(NL_DS_Table) :: table
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_86')
    end if
!
! - Get table
!
    table = ds_energy%table
!
! - Activate columns
!
    do i_col = 1, table%nb_cols
        col_name = table%cols(i_col)%name
        if (col_name .eq. 'NUME_REUSE') then
            table%l_cols_acti(i_col) = ASTER_TRUE
        elseif (col_name .eq. 'INST      ') then
            table%l_cols_acti(i_col) = ASTER_TRUE
        else
!            if (ds_energy%l_comp) then
            table%l_cols_acti(i_col) = ASTER_TRUE
!            end if
        end if
    end do
!
! - Create list of parameters
!
    call nonlinDSTableIOSetPara(table)
!
! - Set other parameters
!
    table%table_io%resultName = resultName
    table%table_io%tablSymbName = 'PARA_CALC'
!
! - Get name of table in results datastructure
!
    call nonlinDSTableIOGetName(table%table_io)
!
! - Create table in results datastructure
!
    call nonlinDSTableIOCreate(table%table_io)
!
! - Save table
!
    ds_energy%table = table
!
end subroutine
