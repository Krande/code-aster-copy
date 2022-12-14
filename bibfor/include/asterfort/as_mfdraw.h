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
    subroutine as_mfdraw(fid, cha, filter, val, locname,&
                         typent, typgeo, numdt, dt, numo,&
                         cret)
        med_idt :: fid
        character(len=*) :: cha
        aster_int :: filter
        real(kind=8) :: val(*)
        character(len=*) :: locname
        aster_int :: typent
        aster_int :: typgeo
        aster_int :: numdt
        real(kind=8) :: dt
        aster_int :: numo
        aster_int :: cret
    end subroutine as_mfdraw
end interface
