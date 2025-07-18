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
subroutine nmextr_read_2(sdextrz, ds_inout, nb_keyw_fact, list_field, rela_field_keyw, &
                         nb_field_comp)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmextd.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
!
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=*), intent(in) :: sdextrz
    integer(kind=8), intent(in) :: nb_keyw_fact
    character(len=24), pointer :: list_field(:)
    integer(kind=8), pointer :: rela_field_keyw(:)
    integer(kind=8), intent(in) :: nb_field_comp
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Read fields to compute
!
! --------------------------------------------------------------------------------------------------
!
! In  sdextr           : name of datastructure for extraction
! In  ds_inout         : datastructure for input/output management
! In  nb_keyw_fact     : number of factor keyword to read extraction parameters
! In  list_field       : list of fields
! In  rela_field_keyw  : relation between field index and keyword index
! In  nb_field_comp    : number of fields to compute (not a default in nonlinear operator)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_keyw_fact, i_field_comp, i_field
    aster_logical :: l_comp
    character(len=2) :: chaine
    character(len=24) :: ligrel
    character(len=24) :: field_type, field_disc, field_comp
    character(len=14) :: sdextr
    character(len=24) :: extr_comp
    character(len=24), pointer :: v_extr_comp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdextr = sdextrz
!
! - Access to datastructure
!
    extr_comp = sdextr(1:14)//'     .COMP'
    call jeveuo(extr_comp, 'E', vk24=v_extr_comp)
!
    do i_field_comp = 1, nb_field_comp
!
! ----- Name of LIGREL
!
        write (chaine, '(I2)') i_field_comp
        ligrel = sdextr(1:14)//chaine(1:2)//'   .LIGR'
!
! ----- Find first keyword for this field
!
        do i_keyw_fact = 1, nb_keyw_fact
            i_field = rela_field_keyw(i_keyw_fact)
            l_comp = i_field .lt. 0
            i_field = abs(i_field)
            if (l_comp) then
                field_type = list_field(i_field)
                if (field_type .eq. 'EPSI_ELGA') then
                    field_disc = 'ELGA'
                    call nmextd(field_type, ds_inout, field_comp)
                else
                    ASSERT(.false.)
                end if
            end if
        end do
!
! ----- Save
!
        v_extr_comp(4*(i_field_comp-1)+1) = field_comp
        v_extr_comp(4*(i_field_comp-1)+2) = field_disc
        v_extr_comp(4*(i_field_comp-1)+3) = field_type
        v_extr_comp(4*(i_field_comp-1)+4) = ligrel
    end do
!
end subroutine
