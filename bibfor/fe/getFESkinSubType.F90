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
subroutine getFESkinSubType(typema, side, subtype, node_init)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrfno.h"
!
    character(len=8), intent(in) :: typema
    character(len=8), intent(in) :: side
    character(len=8), intent(out) :: subtype
    integer(kind=8), intent(out) :: node_init
!
! --------------------------------------------------------------------------------------------------
!
! FE_topo
!
! Get sub type for composed cells like for contact and coupling
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: mast, slav
!
    select case (typema)
    case ("S22")
        slav = "SE2"
        mast = "SE2"
    case ("S32")
        slav = "SE3"
        mast = "SE2"
    case ("S23")
        slav = "SE2"
        mast = "SE3"
    case ("S33")
        slav = "SE3"
        mast = "SE3"
    case ("T33")
        slav = "TR3"
        mast = "TR3"
    case ("T66")
        slav = "TR6"
        mast = "TR6"
    case ("T77")
        slav = "TR7"
        mast = "TR7"
    case ("Q44")
        slav = "QU4"
        mast = "QU4"
    case ("Q88")
        slav = "QU8"
        mast = "QU8"
    case ("Q99")
        slav = "QU9"
        mast = "QU9"
    case ("TQ7")
        slav = "TR3"
        mast = "QU4"
    case ("QT7")
        slav = "QU4"
        mast = "TR3"
    case ("TT1")
        slav = "TR6"
        mast = "TR3"
    case ("TT2")
        slav = "TR3"
        mast = "TR6"
    case ("TQ1")
        slav = "TR6"
        mast = "QU4"
    case ("QT1")
        slav = "QU4"
        mast = "TR6"
    case ("TQ2")
        slav = "TR6"
        mast = "QU8"
    case ("QT2")
        slav = "QU8"
        mast = "TR6"
    case ("TQ3")
        slav = "TR6"
        mast = "QU9"
    case ("QT3")
        slav = "QU9"
        mast = "TR6"
    case ("TQ4")
        slav = "TR3"
        mast = "QU8"
    case ("QT4")
        slav = "QU8"
        mast = "TR3"
    case ("QQ1")
        slav = "QU8"
        mast = "QU4"
    case ("QQ2")
        slav = "QU4"
        mast = "QU8"
    case ("QQ3")
        slav = "QU8"
        mast = "QU9"
    case ("QQ4")
        slav = "QU9"
        mast = "QU8"
    case ("QQ5")
        slav = "QU9"
        mast = "QU4"
    case ("QQ6")
        slav = "QU4"
        mast = "QU9"
    case ("QT5")
        slav = "QU9"
        mast = "TR3"
    case ("TQ5")
        slav = "TR3"
        mast = "QU9"
    case ("QQ7")
        slav = "QU4"
        mast = "QU4"
    case ("QQ8")
        slav = "QU8"
        mast = "QU8"
    case ("TT5")
        slav = "TR6"
        mast = "TR6"
    case ("QT8")
        slav = "QU4"
        mast = "TR7"
    case ("QT9")
        slav = "QU8"
        mast = "TR7"
    case ("QT0")
        slav = "QU9"
        mast = "TR7"
    case ("TT6")
        slav = "TR3"
        mast = "TR7"
    case ("TT7")
        slav = "TR6"
        mast = "TR7"
    case default
        ASSERT(ASTER_FALSE)
    end select
!
    if (side(1:5) == "SLAVE") then
        subtype = slav
        node_init = 1
    elseif (side(1:6) == "MASTER") then
        subtype = mast
        call elrfno(slav, nno=node_init)
        node_init = node_init+1
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
