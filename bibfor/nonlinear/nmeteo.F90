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

subroutine nmeteo(result, sddisc, ds_inout, force, nume_store, &
                  time, i_field, ds_print_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/diincl.h"
#include "asterfort/exisd.h"
#include "asterfort/nmarcc.h"
#include "asterfort/nmetnc.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: sddisc
    character(len=8), intent(in) :: result
    integer(kind=8), intent(in) :: i_field
    integer(kind=8), intent(in) :: nume_store
    real(kind=8), intent(in) :: time
    aster_logical, intent(in) :: force
    type(NL_DS_Print), optional, intent(in) :: ds_print_
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Save field in results datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of datastructure for results
! In  ds_inout         : datastructure for input/output management
! In  nume_store       : index to store in results
! In  i_field          : field index
! In  ds_print         : datastructure for printing parameters
! In  sddisc           : datastructure for discretization
! In  time             : current time
! In  force            : .true. to store field whatever storing options
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=24) :: algo_name, field_algo, field_type
    aster_logical :: l_print, l_store, l_acti
!
! --------------------------------------------------------------------------------------------------
!
!
! - Field to store ?
!
    l_store = ds_inout%field(i_field)%l_store
    if (l_store) then
!
! ----- Print for this step ?
!
        l_print = .true.
        if (present(ds_print_)) then
            l_print = ds_print_%l_print
        end if
!
! ----- Is field should been active ?
!
        l_acti = ds_inout%l_field_acti(i_field)
!
! ----- Name of field (type) in results datastructure
!
        field_type = ds_inout%field(i_field)%type
!
! ----- Name of field in algorithm
!
        algo_name = ds_inout%field(i_field)%algo_name
        call nmetnc(algo_name, field_algo)
!
! ----- Store field
!
        if (l_acti) then
            call exisd('CHAMP', field_algo, iret)
            if (diincl(sddisc, field_type, force) .and. (iret .eq. 1)) then
                if (l_print) then
                    call utmess('I', 'ARCHIVAGE_6', sk=field_type, si=nume_store, sr=time)
                end if
                call nmarcc(result, nume_store, field_type, field_algo)
            end if
        end if
    end if
!
end subroutine
