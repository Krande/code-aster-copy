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
subroutine nmetcc(field_type, algo_name, init_name, &
                  compor, sddyna, ds_contact, &
                  hydr, temp_init, hydr_init)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/ndynkk.h"
!
    character(len=24), intent(in) :: field_type
    character(len=24), intent(out) :: algo_name
    character(len=24), intent(out) :: init_name
    type(NL_DS_Contact), optional, intent(in) :: ds_contact
    character(len=19), optional, intent(in) :: compor
    character(len=19), optional, intent(in) :: sddyna
    character(len=24), optional, intent(in) :: hydr
    character(len=24), optional, intent(in) :: hydr_init
    character(len=24), optional, intent(in) :: temp_init
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Get name of field during non-linear algorithm and initial state
!
! --------------------------------------------------------------------------------------------------
!
! In  field_type       : type of field (symbolic name in result datastructure)
! Out algo_name        : name of field during non-linear algorithm
!                       If 'XXXXXXXXXXXXXXXX' => already defined during DS creation
! Out init_name        : name of field for initial state
!                       If 'XXXXXXXXXXXXXXXX' => already defined during DS creation
! In  compor           : name of <CARTE> COMPOR
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of dynamic parameters datastructure
! In  hydr             : name of field for hydration (HYDR_ELNO)
! In  hydr_init        : name of field for initial hydration
! In  temp_init        : name of field for initial temperature
!
!     if algo_name = #H#TYPCHA#
!       => field name is the TYPCHA "hat" variable datastructure
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: depabs, vitabs, accabs
!
! --------------------------------------------------------------------------------------------------
!
    if (present(sddyna)) then
        call ndynkk(sddyna, 'DEPABS', depabs)
        call ndynkk(sddyna, 'VITABS', vitabs)
        call ndynkk(sddyna, 'ACCABS', accabs)
    end if

! - Standard fields
    if (field_type .eq. 'COMPORTEMENT') then
        algo_name = compor
        init_name = ' '
    else if (field_type .eq. 'CONT_NOEU') then
        algo_name = ds_contact%field_cont_node
        init_name = ' '
    else if (field_type .eq. 'CONT_ELEM') then
        algo_name = ds_contact%field_cont_elem
        init_name = ' '
    else if (field_type .eq. 'DEPL_ABSOLU') then
        algo_name = depabs
        init_name = '&&CNPART.ZERO'
    else if (field_type .eq. 'VITE_ABSOLU') then
        algo_name = vitabs
        init_name = '&&CNPART.ZERO'
    else if (field_type .eq. 'ACCE_ABSOLU') then
        algo_name = accabs
        init_name = '&&CNPART.ZERO'
    else if (field_type .eq. 'TEMP') then
        algo_name = 'XXXXXXXXXXXXXXXX'
        init_name = temp_init
    else if (field_type .eq. 'HYDR_ELGA') then
        algo_name = hydr
        init_name = hydr_init
    else if (field_type .eq. 'COMPORTHER') then
        algo_name = compor
        init_name = ' '
    else
        algo_name = 'XXXXXXXXXXXXXXXX'
        init_name = 'XXXXXXXXXXXXXXXX'
    end if
!
end subroutine
