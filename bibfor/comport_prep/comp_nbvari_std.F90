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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_nbvari_std(rela_comp, defo_comp, type_cpla, &
                           kit_comp, post_iter, regu_visc, post_incr, &
                           nbVari, numeLaw)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
#include "asterfort/comp_meca_code.h"
!
    character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla
    character(len=16), intent(in) :: kit_comp(4), post_iter, regu_visc, post_incr
    integer(kind=8), intent(out) :: nbVari, numeLaw
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Get number of internal state variables for standard constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  type_cpla        : plane stress method
! In  kit_comp         : KIT comportment
! In  post_iter        : type of post-treatment at each Newton iteration
! In  regu_visc        : keyword for viscuous regularization
! In  post_incr        : type of post-treatment at end of time step
! Out nbVari           : number of internal state variables
! Out numeLaw          : index of subroutine for behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: comp_code_py
    integer(kind=8) :: idummy
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = 0
    numeLaw = 0

! - Coding composite comportment (Python)
    call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp, &
                        post_iter, regu_visc, post_incr, &
                        comp_code_py)

! - Get number of total internal state variables and index of law
    call lcinfo(comp_code_py, numeLaw, nbVari, idummy)

! - End of encoding
    call lcdiscard(comp_code_py)
!
end subroutine
