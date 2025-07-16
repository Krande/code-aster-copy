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
subroutine ntetcr(nume_dof, ds_inout, &
                  listLoad_, compor_, hydr_, hydr_init_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nthydr.h"
#include "asterfort/nmetcc.h"
#include "asterfort/vtcreb.h"
#include "asterfort/copisd.h"
#include "asterfort/SetIOField.h"
!
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_InOut), intent(inout) :: ds_inout
    character(len=24), optional, intent(in) :: listLoad_
    character(len=*), optional, intent(in) :: compor_
    character(len=*), optional, intent(in) :: hydr_
    character(len=*), optional, intent(in) :: hydr_init_
!
! --------------------------------------------------------------------------------------------------
!
! THER_* - Init
!
! Create input/output datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_dof         : name of nume_dof object (numbering equation)
! In  compor           : name of <CARTE> COMPOR
! In  listLoad         : name of datastructure for list of loads
! In  hydr             : name of field for hydration
! In  hydr_init        : name of field for initialhydration
! IO  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_field, i_field
    aster_logical :: l_hydr, l_temp_nonl
    character(len=24) :: temp_init, listLoadResu
    character(len=24) :: field_type, algo_name, init_name
    character(len=19) :: compor
    character(len=24) :: hydr, hydr_init
!
! --------------------------------------------------------------------------------------------------
!
    compor = ' '
    hydr = ' '
    hydr_init = ' '
    if (present(compor_)) then
        compor = compor_
    end if
    if (present(hydr_)) then
        hydr = hydr_
    end if
    if (present(hydr_init_)) then
        hydr_init = hydr_init_
    end if
    nb_field = ds_inout%nb_field
    listLoadResu = ds_inout%listLoadResu
    temp_init = '&&NTETCR.TEMP0'
    l_temp_nonl = ds_inout%l_temp_nonl
!
! - Copy of list of loads for save in results datastructure
!
    if (present(listLoad_)) then
        call copisd('LISTE_CHARGES', 'G', listLoad_, listLoadResu)
    end if
!
! - Active functionnalities
!
    l_hydr = .false.
    if (l_temp_nonl) then
        call nthydr(l_hydr)
    end if
!
! - Select fields depending on active functionnalities
!
    call SetIOField(ds_inout, 'TEMP', l_acti_=.true._1)
    if (l_temp_nonl) then
        call SetIOField(ds_inout, 'COMPORTHER', l_acti_=.true._1)
    end if
    if (l_hydr) then
        call SetIOField(ds_inout, 'HYDR_ELGA', l_acti_=.true._1)
    end if
!
! - Add fields
!
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        call nmetcc(field_type, algo_name, init_name, &
                    compor=compor, &
                    hydr=hydr, temp_init=temp_init, hydr_init=hydr_init)
        if (algo_name .ne. 'XXXXXXXXXXXXXXXX') then
            ds_inout%field(i_field)%algo_name = algo_name
            ds_inout%field(i_field)%init_name = init_name
        end if
    end do
!
! - Create initial state fields
!
    call vtcreb(temp_init, 'V', 'R', nume_ddlz=nume_dof)
!
end subroutine
