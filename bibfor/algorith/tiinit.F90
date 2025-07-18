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

subroutine tiinit(ds_inout, sddisc, lostat, l_evol)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ntcrar.h"
#include "asterfort/ntcrli.h"
#include "asterfort/getvid.h"
#include "asterfort/ntcra0.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: sddisc
    aster_logical, intent(in) :: lostat
    aster_logical, intent(out) :: l_evol
!
! --------------------------------------------------------------------------------------------------
!
! THER_* - Datastructures
!
! Time discretization and storing datastructures
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_inout         : datastructure for input/output management
! In  sddisc           : datastructure for time discretization
! In  lostat           : .true. for initial stationnary computation
! Out l_evol           : .true. for transient computation
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: result
    character(len=19) :: listInst
    integer(kind=8) :: nocc
    aster_logical :: l_reuse
!
! --------------------------------------------------------------------------------------------------
!
    l_evol = .false._1
    result = ds_inout%result
    l_reuse = ds_inout%l_reuse
!
! - Transient computation ?
!
    call getvid('INCREMENT', 'LIST_INST', iocc=1, scal=listInst, nbret=nocc)
    if (nocc .eq. 0) then
        if (.not. lostat) then
            call utmess('F', 'DISCRETISATION_8')
        end if
        l_evol = .false.
    else
        l_evol = .true.
    end if

! - Create time discretization datastructure and storing datastructure
    if (l_evol) then
        call ntcrli(listInst, sddisc, lostat)
        call ntcrar(result, sddisc, l_reuse)
    else
        call ntcra0(sddisc)
    end if
!
end subroutine
