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
#include "asterf_types.h"
!
interface
    subroutine xfem_count_no(neq, deeq, k8cmp, nbnomax, ino_xfem, is_xfem, nbnoxfem)
        character(len=8) :: k8cmp(*)
        integer(kind=8) :: deeq(*)
        integer(kind=8) :: neq
        integer(kind=8) :: ino_xfem(nbnomax)
        aster_logical :: is_xfem(nbnomax)
        integer(kind=8) :: nbnomax
        integer(kind=8) :: nbnoxfem
    end subroutine xfem_count_no
end interface
