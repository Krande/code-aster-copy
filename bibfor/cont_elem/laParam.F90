! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine laParam(proj_tole, gamma_c_nodes)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "jeveux.h"
!
real(kind=8), intent(out) :: proj_tole, gamma_c_nodes(4)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get paramaters from mmchml_a.F90
!
! --------------------------------------------------------------------------------------------------
!
    integer :: jcont
!
    call jevech('PCONFR', 'L', jcont)
!
! - Parameters
!
    proj_tole = zr(jcont-1+23)
    gamma_c_nodes = zr(jcont-1+24)
!
end subroutine
