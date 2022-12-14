! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine ntetcr(nume_dof  , ds_inout,&
                  list_load_, compor_ , hydr_, hydr_init_)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nthydr.h"
#include "asterfort/nmetcc.h"
#include "asterfort/vtcreb.h"
#include "asterfort/copisd.h"
#include "asterfort/SetIOField.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_InOut), intent(inout) :: ds_inout
    character(len=19), optional, intent(in) :: list_load_
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
! In  list_load        : name of datastructure for list of loads
! In  hydr             : name of field for hydratation
! In  hydr_init        : name of field for initialhydratation
! IO  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_field, i_field
    aster_logical :: l_hydr, l_temp_nonl
    character(len=24) :: temp_init
    character(len=24) :: field_type, algo_name, init_name
    character(len=19) :: compor, list_load_resu
    character(len=24) :: hydr, hydr_init
!
! --------------------------------------------------------------------------------------------------
!
    compor       = ' '
    hydr         = ' '
    hydr_init    = ' '
    if (present(compor_))    then
        compor = compor_
    endif
    if (present(hydr_))      then
        hydr = hydr_
    endif
    if (present(hydr_init_)) then
        hydr_init = hydr_init_
    endif
    nb_field       = ds_inout%nb_field
    list_load_resu = ds_inout%list_load_resu
    temp_init      = '&&NTETCR.TEMP0'
    l_temp_nonl    = ds_inout%l_temp_nonl
!
! - Copy of list of loads for save in results datastructure
!
    if (present(list_load_)) then
        call copisd('LISTE_CHARGES', 'G', list_load_, list_load_resu)
    endif
!
! - Active functionnalities
!
    l_hydr = .false.
    if (l_temp_nonl) then
        call nthydr(l_hydr)
    endif
!
! - Select fields depending on active functionnalities
!
    call SetIOField(ds_inout, 'TEMP', l_acti_ = .true._1)
    if (l_temp_nonl) then
        call SetIOField(ds_inout, 'COMPORTHER', l_acti_ = .true._1)
    endif
    if (l_hydr) then
        call SetIOField(ds_inout, 'HYDR_ELNO' , l_acti_ = .true._1)
    endif
!
! - Add fields
!
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        call nmetcc(field_type     , algo_name, init_name, &
                    compor = compor ,&
                    hydr   = hydr    , temp_init = temp_init, hydr_init = hydr_init)
        if (algo_name.ne.'XXXXXXXXXXXXXXXX') then
            ds_inout%field(i_field)%algo_name = algo_name
            ds_inout%field(i_field)%init_name = init_name
        endif
    end do
!
! - Create initial state fields
!
    call vtcreb(temp_init, 'V', 'R', nume_ddlz = nume_dof)
!
end subroutine
