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
subroutine nmimen(ds_print)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nonlinDSColumnWriteValue.h"
#include "asterfort/iunifi.h"
#include "asterfort/ComputeTableHead.h"
!
    type(NL_DS_Print), intent(in) :: ds_print
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Print management
!
! Print head of convergence table
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_print         : datastructure for printing parameters
!
! --------------------------------------------------------------------------------------------------
!
    character(len=512) :: sep_line, table_head(3)
    type(NL_DS_Table) :: table_cvg
    aster_logical :: l_tcvg_csv
    integer(kind=8) :: line_width, mesg_unit, tcvg_unit
!
! --------------------------------------------------------------------------------------------------
!
    mesg_unit = iunifi('MESSAGE')
!
! - Get convergence table
!
    table_cvg = ds_print%table_cvg
!
! - Get parameters
!
    sep_line = ds_print%table_cvg%sep_line
    tcvg_unit = ds_print%tcvg_unit
    line_width = ds_print%table_cvg%width
    l_tcvg_csv = ds_print%l_tcvg_csv
!
! - Compute table head
!
    call ComputeTableHead(table_cvg, '|', table_head)
!
! - Print in message unit
!
    call nonlinDSColumnWriteValue(line_width, &
                                  output_unit_=mesg_unit, &
                                  value_k_=sep_line)
    call nonlinDSColumnWriteValue(line_width, &
                                  output_unit_=mesg_unit, &
                                  value_k_=table_head(1))
    call nonlinDSColumnWriteValue(line_width, &
                                  output_unit_=mesg_unit, &
                                  value_k_=table_head(2))
    call nonlinDSColumnWriteValue(line_width, &
                                  output_unit_=mesg_unit, &
                                  value_k_=table_head(3))
    call nonlinDSColumnWriteValue(line_width, &
                                  output_unit_=mesg_unit, &
                                  value_k_=sep_line)
!
! - Print in file
!
    if (l_tcvg_csv) then
        call ComputeTableHead(table_cvg, ',', table_head)
        call nonlinDSColumnWriteValue(line_width, &
                                      output_unit_=tcvg_unit, &
                                      value_k_=table_head(1))
        call nonlinDSColumnWriteValue(line_width, &
                                      output_unit_=tcvg_unit, &
                                      value_k_=table_head(2))
        call nonlinDSColumnWriteValue(line_width, &
                                      output_unit_=tcvg_unit, &
                                      value_k_=table_head(3))
    end if
!
end subroutine
