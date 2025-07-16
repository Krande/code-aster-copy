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

!
!
#include "asterf_types.h"
interface
    subroutine asmpi_sendrecv_i(buffer_send, count_send, recipient, tag_send, &
                                buffer_recv, count_recv, sender   , tag_recv, &
                                comm)
        integer(kind=8), intent(in) :: buffer_send(*)
        integer(kind=8), intent(out) :: buffer_recv(*)
        mpi_int, intent(in) :: count_send, count_recv
        mpi_int, intent(in) :: recipient, sender
        mpi_int, intent(in) :: tag_send, tag_recv
        mpi_int, intent(in) :: comm
    end subroutine asmpi_sendrecv_i
end interface
