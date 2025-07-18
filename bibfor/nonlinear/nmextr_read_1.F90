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

subroutine nmextr_read_1(ds_inout, keyw_fact, nb_keyw_fact, list_field, rela_field_keyw, &
                         nb_field, nb_field_comp)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/nmextc.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_InOut), intent(in) :: ds_inout
    integer(kind=8), intent(in) :: nb_keyw_fact
    character(len=16), intent(in) :: keyw_fact
    character(len=24), pointer :: list_field(:)
    integer(kind=8), pointer :: rela_field_keyw(:)
    integer(kind=8), intent(out) :: nb_field
    integer(kind=8), intent(out) :: nb_field_comp
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Read fields to extract
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_inout         : datastructure for input/output management
! In  keyw_fact        : factor keyword to read extraction parameters
! In  nb_keyw_fact     : number of factor keyword to read extraction parameters
! Out list_field       : list of fields
! Out rela_field_keyw  : relation between field index and keyword index
! Out nb_field         : total number of fields
! Out nb_field_comp    : number of fields to compute (not a default in nonlinear operator)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_keyw_fact, i_field, i_list_field
    character(len=24) :: field_type, field_old
    aster_logical :: l_extr, l_find
!
! --------------------------------------------------------------------------------------------------
!
    nb_field = 0
    nb_field_comp = 0
    if (nb_keyw_fact .eq. 0) then
        goto 99
    end if
!
! - List of field to extract
!
    AS_ALLOCATE(vk24=list_field, size=nb_keyw_fact)
!
! - Relation between field index and keyword index
!
    AS_ALLOCATE(vi=rela_field_keyw, size=nb_keyw_fact)
!
    do i_keyw_fact = 1, nb_keyw_fact
!
! ----- Read field type
!
        call nmextc(ds_inout, keyw_fact, i_keyw_fact, field_type, l_extr)
        if (.not. l_extr) then
            field_type = 'NONE'
        end if
!
! ----- Add field in list to extract
!
        l_find = .false.
        do i_list_field = 1, nb_keyw_fact-1
            field_old = list_field(i_list_field)
            if (field_old .eq. field_type) then
                i_field = i_list_field
                l_find = .true.
            end if
        end do
        if (.not. l_find) then
            nb_field = nb_field+1
            i_field = nb_field
            list_field(i_field) = field_type
            if (field_type .eq. 'EPSI_ELGA') then
                nb_field_comp = nb_field_comp+1
            end if
        end if
!
! ----- Set relation between field index and keyword index
!
        if (field_type .eq. 'EPSI_ELGA') then
            rela_field_keyw(i_keyw_fact) = -i_field
        else
            rela_field_keyw(i_keyw_fact) = i_field
        end if
    end do
!
99  continue
!
end subroutine
