! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine evenPrtIncrQuan(mesgOrigZ, sddisc, index, newdt_)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: mesgOrigZ
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: index
    real(kind=8), optional, intent(in) :: newdt_
!
! --------------------------------------------------------------------------------------------------
!
! Management of event
!
! Print and save informations about DELTA_GRANDEUR
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4) :: mesgOrig
    integer(kind=8) ::  iAdap, iEven
    character(len=16) :: fieldType, cmpName
!
! --------------------------------------------------------------------------------------------------
!
    mesgOrig = mesgOrigZ
    if (mesgOrig .eq. "ERRE") then
        iEven = index
        call utdidt('L', sddisc, 'ECHE', 'NOM_CHAM', index_=iEven, valk_=fieldType)
        call utdidt('L', sddisc, 'ECHE', 'NOM_CMP', index_=iEven, valk_=cmpName)
        call utmess('I', 'MECANONLINE10_24', nk=2, valk=[fieldType, cmpName])

    elseif (mesgOrig .eq. "ADAP") then
        iAdap = index
        call utdidt('L', sddisc, 'ADAP', 'NOM_CHAM', index_=iAdap, valk_=fieldType)
        call utdidt('L', sddisc, 'ADAP', 'NOM_CMP', index_=iAdap, valk_=cmpName)
        call utmess('I', 'ADAPTATION_20', nk=2, valk=[fieldType, cmpName], sr=newdt_)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
