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
subroutine pmevdr(sddisc, tablType, tablIncr, &
                  liccvg, lIterNewtMaxi, conver, &
                  newtLoopAction)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/nmacto.h"
#include "asterfort/pmevel.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: tablType
    character(len=24), intent(in) :: tablIncr
    integer(kind=8), intent(in) :: liccvg(5)
    aster_logical, intent(in) :: lIterNewtMaxi, conver
    integer(kind=8), intent(out) :: newtLoopAction
!
! --------------------------------------------------------------------------------------------------
!
! SIMU_POINT_MAT
!
! Detect event
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: faccvg, ldccvg
    aster_logical :: lerror
    integer(kind=8) :: ievdac
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<SIMUPOINTMAT> EVALUATION DES EVENT-DRIVEN'
    end if

! - INITIALISATIONS
    ldccvg = liccvg(2)
    faccvg = liccvg(5)
    lerror = (ldccvg .eq. 1) .or. (faccvg .ne. 0) .or. lIterNewtMaxi

! - NEWTON A CONVERGE ?
    if (conver) then
        newtLoopAction = 0
    else
        newtLoopAction = 2
    end if

! - Detect first event
    call pmevel(sddisc, tablType, tablIncr, &
                lerror, conver)

! - UN EVENEMENT SE DECLENCHE
    call nmacto(sddisc, ievdac)
    if (ievdac .ne. 0) then
        newtLoopAction = 1
    end if
!
end subroutine
