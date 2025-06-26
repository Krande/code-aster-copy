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
! person_in_charge: nicolas.pignet at edf.fr
!
function asmpi_all(ibool, test) result(obool)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
!
    aster_logical, intent(in) :: ibool, test
    aster_logical             :: obool
!
!---------------------------------------------------------------------------------------------------
!   But :
!     To know if the boolean is true or false on all processes
!
!   IN:
!     ibool   : local input boolean
!     test    : value to check
!
!   OUT:
!     obool   : gloab output boolean
!
!---------------------------------------------------------------------------------------------------
    integer(kind=8) :: iint
!-----------------------------------------------------------------------
!
    if (ibool .eqv. test) then
        iint = 0
    else
        iint = 1
    end if
!
    call asmpi_comm_vect('MPI_MAX', 'I', sci=iint)

    if (iint == 0) then
        obool = ASTER_TRUE
    else
        obool = ASTER_FALSE
    end if
!
end function
