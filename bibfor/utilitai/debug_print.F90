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

subroutine debug_print(sch1, unit)
    implicit none
!     --
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/utimsd.h"
    character(len=*), intent(in) :: sch1
    integer(kind=8), optional, intent(in) :: unit
!
    integer(kind=8) :: rang, nbproc, iproc, unit_
    mpi_int :: mrank, msize
!
    if (present(unit)) then
        unit_ = unit
    else
        unit_ = 6
    end if
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    do iproc = 0, nbproc-1
        if (iproc .eq. rang) then
            call utimsd(unit_, 2, .false._1, .true._1, sch1, 1, ' ')
            flush (unit_)
        end if
        call asmpi_barrier()
    end do
end subroutine
