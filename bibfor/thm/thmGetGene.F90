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
subroutine thmGetGene(ds_thm, l_vf, ndim, &
                      mecani, press1, press2, tempe, second)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(THM_DS), intent(in) :: ds_thm
    aster_logical, intent(in) :: l_vf
    integer(kind=8), intent(in)  :: ndim
    integer(kind=8), intent(out) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
!
! --------------------------------------------------------------------------------------------------
!
! THM - Initializations
!
! Get generalized coordinates
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  l_vf             : flag for finite volume
! In  ndim             : dimension of space (2 or 3)
! Out mecani           : parameters for mechanic
!                    (1) - Flag if physic exists (1 if exists)
!                    (2) - Adress of first component in generalized strain vector
!                    (3) - Adress of first component in generalized stress vector
!                    (4) - Number of components for strains
!                    (5) - Number of components for stresses
! Out press1           : parameters for hydraulic (first pressure)
!                    (1) - Flag if physic exists (1 if exists)
!                    (2) - Number of phases
!                    (3) - Adress of first component in generalized strain vector
!                    (4) - Adress of first component in vector of gen. stress for first phase
!                    (5) - Adress of first component in vector of gen. stress for second phase
!                    (6) - Number of components for strains
!                    (7) - Number of components for stresses (for each phase)
! Out press1           : parameters for hydraulic (second pressure)
!                    (1) - Flag if physic exists (1 if exists)
!                    (2) - Number of phases
!                    (3) - Adress of first component in generalized strain vector
!                    (4) - Adress of first component in vector of gen. stress for first phase
!                    (5) - Adress of first component in vector of gen. stress for second phase
!                    (6) - Number of components for strains
!                    (7) - Number of components for stresses (for each phase)
! Out tempe            : parameters for thermic
!                    (1) - Flag if physic exists (1 if exists)
!                    (2) - Adress of first component in generalized strain vector
!                    (3) - Adress of first component in generalized stress vector
!                    (4) - Number of components for strains
!                    (5) - Number of components for stresses
! Out second           : parameters for second gradient
!                    (1) - Flag if physic exists (1 if exists)
!                    (2) - Adress of first component in generalized strain vector
!                    (3) - Adress of first component in generalized stress vector
!                    (4) - Number of components for strains
!                    (5) - Number of components for stresses
! --------------------------------------------------------------------------------------------------
!
    mecani(:) = 0
    press1(:) = 0
    press2(:) = 0
    tempe(:) = 0
    second(:) = 0
!
! - Main parameters: mechanic, thermic, hydraulic
!
    if (ds_thm%ds_elem%l_dof_meca) then
        mecani(1) = 1
    end if
    if (ds_thm%ds_elem%l_dof_ther) then
        tempe(1) = 1
    end if
    if (ds_thm%ds_elem%l_dof_pre1) then
        press1(1) = 1
    end if
    if (ds_thm%ds_elem%l_dof_pre2) then
        press2(1) = 1
    end if
    if (ds_thm%ds_elem%l_dof_2nd) then
        second(1) = 1
    end if
    press1(2) = ds_thm%ds_elem%nb_phase(1)
    press2(2) = ds_thm%ds_elem%nb_phase(2)
!
! - Number of (generalized) stress/strain components - Mechanic
!
    if (mecani(1) .eq. 1) then
        mecani(4) = ndim+6
        mecani(5) = 6+6
    end if
!
! - Number of (generalized) stress/strain components - Thermic
!
    if (tempe(1) .eq. 1) then
        tempe(4) = 1+ndim
        tempe(5) = 1+ndim
    end if
!
! - Number of (generalized) stress/strain components - Hydraulic
!
    if (press1(1) .eq. 1) then
        press1(6) = 1+ndim
        press1(7) = ndim+1
        if (tempe(1) .eq. 1) then
            press1(7) = press1(7)+1
        end if
    end if
    if (press2(1) .eq. 1) then
        press2(6) = 1+ndim
        press2(7) = ndim+1
        if (tempe(1) .eq. 1) then
            press2(7) = press2(7)+1
        end if
    end if
    if (l_vf) then
        press1(7) = 2
        press2(7) = 2
    end if
!
! - Number of (generalized) stress/strain components - Second gradient
!
    if (second(1) .eq. 1) then
        second(4) = ndim+3
        second(5) = ndim+3
    end if
!
! - Index for adress in (generalized) vectors - Mechanic
!
    if (mecani(1) .eq. 1) then
        mecani(2) = 1
        mecani(3) = 1
    end if
!
! - Index for adress in (generalized) vectors - Hydraulic
!
    if (press1(1) .eq. 1) then
        press1(3) = mecani(4)+1
        press1(4) = mecani(5)+1
        if (press1(2) .eq. 2) then
            press1(5) = press1(4)+press1(7)
        end if
    end if
!
    if (press2(1) .eq. 1) then
        press2(3) = press1(3)+press1(6)
        press2(4) = press1(4)+press1(2)*press1(7)
        if (press2(2) .eq. 2) then
            press2(5) = press2(4)+press2(7)
        end if
    end if
!
! - Index for adress in (generalized) vectors - Thermic
!
    if (tempe(1) .eq. 1) then
        tempe(2) = mecani(4)+press1(6)+press2(6)+1
        tempe(3) = mecani(5)+press1(2)*press1(7)+press2(2)*press2(7)+1
    end if
!
! - Index for adress in (generalized) vectors - Second gradient
!
    if (second(1) .eq. 1) then
        second(2) = mecani(4)+press1(6)+press2(6)+tempe(4)+1
        second(3) = mecani(5)+press1(2)*press1(7)+press2(2)*press2(7)+tempe(5)+1
    end if
!

end subroutine
