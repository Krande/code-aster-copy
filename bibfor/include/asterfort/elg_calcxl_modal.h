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
#include "asterf_petsc.h"
!
interface
    subroutine elg_calcxl_modal(x1, omega2, ke_mass, vlag)
#ifdef ASTER_HAVE_PETSC
        use aster_petsc_module
        Vec :: x1, vlag
        PetscInt :: ke_mass
        real(kind=8) :: omega2
# else
        integer(kind=8) :: x1, vlag, ke_mass
        real(kind=8) :: omega2
# endif
    end subroutine elg_calcxl_modal
end interface
