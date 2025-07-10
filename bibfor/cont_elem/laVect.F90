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
subroutine laVect(parameters, geom, vect_cont, vect_fric)
!
    use contact_module
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/laVect_ct_pr.h"
#include "asterfort/laVect_ct_std.h"
#include "asterfort/laVect_cf_pr.h"
#include "asterfort/laVect_cf_std.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: vect_cont(MAX_LAGA_DOFS), vect_fric(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute vector
!
! -----
!
! Switch between different case : frictionless contact, frictional contact following
! Poulios & Renard's formulation, frictional contact following classical formulation
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : numerical parameters
! In  geom             : geometrical information
! IO  vect_cont        : vector for contact
! IO  vect_fric        : vector for friction
!
! --------------------------------------------------------------------------------------------------
!
!
    if (.not. parameters%l_fric) then
        if (parameters%vari_cont == CONT_VARI_ROBU .or. &
            parameters%vari_cont == CONT_VARI_NONE) then
            ! Frictionless contact following Poulios & Renard's formulation
            call laVect_ct_pr(parameters, geom, vect_cont, vect_fric)
            ! write (6, *) 'laVect_ct_pr'
!
        elseif (parameters%vari_cont == CONT_VARI_CLAS) then
            ! Frictionless contact following classical formulation
            call laVect_ct_std(parameters, geom, vect_cont, vect_fric)
            ! write (6, *) 'laVect_ct_std'
        else
            ASSERT(ASTER_FALSE)
        end if
!
    else if (parameters%vari_cont == CONT_VARI_ROBU .or. &
             parameters%vari_cont == CONT_VARI_NONE) then
        ! Frictional contact following Poulios & Renard's formulation
        call laVect_cf_pr(parameters, geom, vect_cont, vect_fric)
        ! write (6, *) 'laVect_cf_pr'
!
    elseif (parameters%vari_cont == CONT_VARI_CLAS) then
        ! Frictional contact following classical formulation
        call laVect_cf_std(parameters, geom, vect_cont, vect_fric)
        ! write (6, *) 'laVect_cf_std'
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
