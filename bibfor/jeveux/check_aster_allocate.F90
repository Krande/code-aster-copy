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

subroutine check_aster_allocate(stage)
    use allocate_module
    implicit none
    integer(kind=8), intent(in) :: stage

! Check that objects allocated by 'as_allocate' have actually been deallocated.
!   stage = 0: Reset the value cuvtrav (from common) to 0 (called before each command)
!   stage = 1: Check and deallocate objects (called after each command)
!   stage = 2: Only dellocate objects (because an exception was raised)
!
#include "jeveux_private.h"
#include "asterc/jdcget.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), save :: icode = -1
!
    if (stage .eq. 0) then
        cuvtrav = 0.d0
    else if (stage .eq. 1) then
        if (abs(cuvtrav) > r8prem()) then
            call utmess('A', 'DVP_6', sr=cuvtrav*lois/1.e6)
        end if
        if (icode < 0) then
            icode = jdcget('TestMode')
        end if
        if (icode .ne. 0) then
            ASSERT(abs(cuvtrav) < r8prem())
        end if
        call deallocate_all_slvec()
    else
        call deallocate_all_slvec()
    end if
!
end subroutine
