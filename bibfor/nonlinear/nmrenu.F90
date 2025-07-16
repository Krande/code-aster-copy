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
subroutine nmrenu(modelz, list_func_acti, list_load, &
                  ds_measure, ds_contact, nume_dof, &
                  l_renumber)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/numer3.h"
#include "asterfort/utmess.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
!
    character(len=*), intent(in) :: modelz
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=19), intent(in) :: list_load
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=24), intent(inout) :: nume_dof
    aster_logical, intent(out) :: l_renumber
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm
!
! Renumbering equations ?
!
! --------------------------------------------------------------------------------------------------
!
! IO  nume_dof         : name of numbering object (NUME_DDL)
! In  model            : name of model datastructure
! In  list_load        : list of loads
! IO  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
! In  list_func_acti   : list of active functionnalities
! Out l_renumber       : .true. if renumber
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=2), parameter :: base = "VG"
    aster_logical :: l_cont, l_cont_cont, l_cont_elem
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)

! - No renumbering !
    l_renumber = ASTER_FALSE

! - Active functionnalities
    l_cont = isfonc(list_func_acti, 'CONTACT')
    l_cont_elem = isfonc(list_func_acti, 'ELT_CONTACT')
    l_cont_cont = isfonc(list_func_acti, 'CONT_CONTINU')

! - To change numbering
    if (l_cont) then
! ----- Start timer for preparation of contact
        call nmtime(ds_measure, 'Launch', 'Cont_Prep')

! ----- Numbering to change ?
        if (l_cont_elem) then
            l_renumber = ds_contact%l_renumber
            ds_contact%l_renumber = ASTER_FALSE
        end if

! ----- Re-numbering
        if (l_renumber) then
            if (niv .ge. 2) then
                call utmess('I', 'MECANONLINE13_36')
            end if
            call numer3(modelZ, base, list_load, nume_dof, ds_contact)
        end if

! ----- Stop timer for preparation of contact
        call nmtime(ds_measure, 'Stop', 'Cont_Prep')
        call nmrinc(ds_measure, 'Cont_Prep')
    end if
!
end subroutine
