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

subroutine cesred_wrap(ces1z, nbma, lima, nbcmp, licmp, base, ces2z)
! person_in_charge: jacques.pellet at edf.fr
! A_UTIL
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cesred.h"
    integer(kind=8) :: nbma, nbcmp
    integer(kind=8) :: lima(nbma)
    character(len=*) :: ces1z, ces2z, base
    character(len=8) :: licmp(*)

    call cesred(ces1z, nbma, lima, nbcmp, licmp, base, ces2z)
end subroutine
