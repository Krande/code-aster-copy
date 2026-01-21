! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
module coupling_type
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/coupling_penalisation_module.h"
#include "asterfort/HHO_size_module.h"
#include "FE_basis_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - Generic
!
! --------------------------------------------------------------------------------------------------
!
!
    type CouplingMap
!
        integer(kind=8) :: nbDoFs = 0, nbDoFsFECEll = 0, nbDoFsFEFace = 0, nbDoFshhoFace = 0
        integer(kind=8) :: mapDoFsFECell(MAX_BV) = 0
        integer(kind=8) :: mapDoFsFEFace(MAX_BV) = 0
        integer(kind=8) :: mapDoFshhoFace(MSIZE_FACE_VEC) = 0
    end type
!
    type CouplingData
!
        real(kind=8) :: disp_prev(MSIZE_FACE_VEC+MAX_BV) = 0.d0
        real(kind=8) :: dispFECell_prev(MAX_BV) = 0.d0
        real(kind=8) :: dispFEFace_prev(MAX_BV) = 0.d0
        real(kind=8) :: disphhoFace_prev(MSIZE_FACE_VEC) = 0.d0
!
        real(kind=8) :: disp_curr(MSIZE_FACE_VEC+MAX_BV) = 0.d0
        real(kind=8) :: dispFECell_curr(MAX_BV) = 0.d0
        real(kind=8) :: dispFEFace_curr(MAX_BV) = 0.d0
        real(kind=8) :: disphhoFace_curr(MSIZE_FACE_VEC) = 0.d0
!
        real(kind=8) :: coef_pena = 0.d0
    end type
!
    public :: CouplingMap, CouplingData
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
end module
