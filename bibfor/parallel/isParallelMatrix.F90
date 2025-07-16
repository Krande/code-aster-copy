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
function isParallelMatrix(matrix) result(l_parallel_matrix)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
!
    character(len=*), intent(in) :: matrix
    aster_logical                :: l_parallel_matrix
!
!---------------------------------------------------------------------------------------------------
!   But :
!     To know if the matrix is parallel_matrix
!
!   IN:
!     matrix      : name of the matrix
!
!   OUT:
!     l_parallel_matrix : the matrix is a parallel_matrix ?
!
!---------------------------------------------------------------------------------------------------
    character(len=3) :: mathpc
    character(len=19) :: matas
!-----------------------------------------------------------------------
!
    matas = matrix
!
    call dismoi('MATR_HPC', matas, 'MATR_ASSE', repk=mathpc)
    l_parallel_matrix = (mathpc == "OUI")
!
end function
