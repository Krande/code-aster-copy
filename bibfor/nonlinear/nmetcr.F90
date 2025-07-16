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
subroutine nmetcr(ds_inout, model, compor, list_func_acti, sddyna, &
                  ds_contact, cara_elem, listLoad)
!
    use NonLin_Datastructure_type
    use listLoad_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/nmetac.h"
#include "asterfort/nmetc0.h"
#include "asterfort/nmetcc.h"
#include "asterfort/rscrsd.h"
!
    type(NL_DS_InOut), intent(inout) :: ds_inout
    character(len=24), intent(in) :: model
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=24), intent(in) :: compor
    character(len=19), intent(in) :: sddyna
    character(len=24), intent(in) :: cara_elem
    character(len=19), intent(in) :: listLoad
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output management
!
! Initializations for input/output management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_inout         : datastructure for input/output management
! In  model            : name of model
! In  cara_elem        : name of datastructure for elementary parameters (CARTE)
! In  compor           : name of <CARTE> COMPOR
! In  ds_contact       : datastructure for contact management
! In  list_func_acti   : list of active functionnalities
! In  sddyna           : name of dynamic parameters datastructure
! In  list_load        : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_field, nb_field_resu
    integer(kind=8) :: i_field, i_field_resu
    aster_logical :: l_find
    character(len=19) :: result
    character(len=24) :: field_resu, field_type, algo_name, init_name, listLoadResu
!
! --------------------------------------------------------------------------------------------------
!
    result = '&&NMETCR'
    nb_field = ds_inout%nb_field
    listLoadResu = ds_inout%listLoadResu

! - Copy of list of loads for save in results datastructure
    call copisd('LISTE_CHARGES', 'G', listLoad, listLoadResu)

! - Select fields depending on active functionnalities
    call nmetac(list_func_acti, sddyna, ds_contact, ds_inout)

! - Add fields
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        call nmetcc(field_type, algo_name, init_name, &
                    compor, sddyna, ds_contact)
        if (algo_name .ne. 'XXXXXXXXXXXXXXXX') then
            ds_inout%field(i_field)%algo_name = algo_name
        end if
        if (init_name .ne. 'XXXXXXXXXXXXXXXX') then
            ds_inout%field(i_field)%init_name = init_name
        end if
    end do
!
! - Create initial state fields
!
    call nmetc0(model, cara_elem, compor, ds_inout)
!
! - Check: fields have been defined in rscrsd.F90 ?
!
    call rscrsd('V', result, 'EVOL_NOLI', 1)
    call jelira(result(1:8)//'           .DESC', 'NOMMAX', nb_field_resu)
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        init_name = ds_inout%field(i_field)%init_name
        if (ds_inout%field(i_field)%l_store) then
            l_find = .false._1
            do i_field_resu = 1, nb_field_resu
                call jenuno(jexnum(result(1:8)//'           .DESC', i_field_resu), field_resu)
                if (field_resu .eq. field_type) then
                    l_find = .true._1
                end if
            end do
! --------- No field in results => change rscrsd subroutine !
            ASSERT(l_find)
        end if
    end do
    call detrsd('RESULTAT', result)
!
end subroutine
