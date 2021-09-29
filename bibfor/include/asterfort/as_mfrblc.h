! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
!
#include "asterf_types.h"
interface
    subroutine as_mfrblc(fid, nent, nvent, ncent, cs,&
                         swm, stm, pname, start, stride,&
                         count, bsize, lbsize, flt, cret)
        med_idt :: fid
        aster_int :: nent
        aster_int :: nvent
        aster_int :: ncent
        aster_int :: cs
        aster_int :: swm
        aster_int :: stm
        aster_int :: start
        aster_int :: stride
        aster_int :: count
        aster_int :: bsize
        aster_int :: lbsize
        aster_int :: flt
        aster_int :: cret
        character(len=*) :: pname
    end subroutine as_mfrblc
end interface
