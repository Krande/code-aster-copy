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

function mminfr(sdcont_defi_, question_, i_zone_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/mminfp.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    real(kind=8) :: mminfr
    character(len=*), intent(in) :: sdcont_defi_
    character(len=*), intent(in) :: question_
    integer(kind=8), optional, intent(in) :: i_zone_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Get parameter (real) - By zone contact
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  question         : question to select parameter
! In  i_zone           : index of contact zone
! Out mminfr           : value for selected parameter
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_zone
    real(kind=8) :: vale_r
    character(len=24) :: sdcont_defi
!
! --------------------------------------------------------------------------------------------------
!
    if (present(i_zone_)) then
        i_zone = i_zone_
    else
        i_zone = 1
    end if
    sdcont_defi = sdcont_defi_
    call mminfp(i_zone, sdcont_defi, question_, vale_r_=vale_r)
    mminfr = vale_r
end function
