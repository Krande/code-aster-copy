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
subroutine nsetcr(nume_dof, ds_inout, &
                  listLoad_, compor_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmetcc.h"
#include "asterfort/vtcreb.h"
#include "asterfort/copisd.h"
#include "asterfort/SetIOField.h"
!
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_InOut), intent(inout) :: ds_inout
    character(len=24), optional, intent(in) :: listLoad_
    character(len=*), optional, intent(in) :: compor_
!
! --------------------------------------------------------------------------------------------------
!
! SECH_NON_LINE - Init
!
! Create input/output datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_dof         : name of nume_dof object (numbering equation)
! In  compor           : name of <CARTE> COMPOR
! In  listLoad         : name of datastructure for list of loads
! IO  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_field, i_field
    aster_logical :: l_hydr, l_temp_nonl
    character(len=24) :: sech_init, listLoadResu
    character(len=24) :: field_type, algo_name, init_name
    character(len=19) :: compor
    character(len=24) :: hydr, hydr_init
!
! --------------------------------------------------------------------------------------------------
!
    compor = ' '
    if (present(compor_)) then
        compor = compor_
    end if
    nb_field = ds_inout%nb_field
    listLoadResu = ds_inout%listLoadResu
    sech_init = '&&NTETCR.SECH0'
    l_temp_nonl = ds_inout%l_temp_nonl
!
! - Copy of list of loads for save in results datastructure
!
    if (present(listLoad_)) then
        call copisd('LISTE_CHARGES', 'G', listLoad_, listLoadResu)
    end if
!
! - Select fields depending on active functionnalities
!
    call SetIOField(ds_inout, 'SECH', l_acti_=.true._1)
    if (l_temp_nonl) then
        call SetIOField(ds_inout, 'COMPORTHER', l_acti_=.true._1)
    end if
!
! - Add fields
!
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        call nmetcc(field_type, algo_name, init_name, &
                    compor=compor, &
                    temp_init=sech_init)
        if (algo_name .ne. 'XXXXXXXXXXXXXXXX') then
            ds_inout%field(i_field)%algo_name = algo_name
            ds_inout%field(i_field)%init_name = init_name
        end if
    end do
!
! - Create initial state fields
!
    call vtcreb(sech_init, 'V', 'R', nume_ddlz=nume_dof)
!
end subroutine
