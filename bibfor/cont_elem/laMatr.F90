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
subroutine laMatr(parameters, geom, matr_cont, matr_fric)
!
    use contact_type
    use contact_module
    use contact_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/laMatr_ct_pr.h"
#include "asterfort/laMatr_ct_std.h"
#include "asterfort/laMatr_cf_pr.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: matr_cont(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    real(kind=8), intent(inout) :: matr_fric(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute matrices
!
! Switch between different case : frictionless contact, frictional contact following
! Poulios & Renard's formulation, frictional contact following classical formulation
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : numerical parameters
! In  geom             : geometrical information
! IO  matr_cont        : matrix (only upper part)
! IO  matr_fric        : matrix
!
! --------------------------------------------------------------------------------------------------
!
    ! Frictionless contact
    if (.not. parameters%l_fric) then
        if (parameters%vari_cont == CONT_VARI_ROBU .or. &
            parameters%vari_cont == CONT_VARI_NONE) then
            ! Frictionless contact following Poulios & Renard's formulation
            call laMatr_ct_pr(parameters, geom, matr_cont, matr_fric)
        elseif (parameters%vari_cont == CONT_VARI_CLAS) then
            ! Frictionless contact following classical formulation
            call laMatr_ct_std(parameters, geom, matr_cont, matr_fric)
        end if
!
    else if (parameters%vari_cont == CONT_VARI_ROBU .or. &
             parameters%vari_cont == CONT_VARI_NONE) then
        ! Frictional contact following Poulios & Renard's formulation
        call laMatr_cf_pr(parameters, geom, matr_cont, matr_fric)
!
    elseif (parameters%vari_cont == CONT_VARI_CLAS) then
        ! Frictional contact following classical formulation
        ! Not implemented yet - finite difference must be used
        ASSERT(ASTER_FALSE)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
