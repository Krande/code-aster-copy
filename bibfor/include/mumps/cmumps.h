! --------------------------------------------------------------------
! Copyright (C) 2007 - 2021 - EDF R&D - www.code-aster.org
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
#ifdef ASTER_HAVE_MUMPS
interface
#ifdef ASTER_PLATFORM_MSVC64
    subroutine cmumps(cmpsk) bind(C, name='CMUMPS')
#else
    subroutine cmumps(cmpsk)
#endif
#       include "cmumps_struc.h"
        type (cmumps_struc) :: cmpsk
    end subroutine cmumps
end interface
#endif
