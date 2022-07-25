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
subroutine laParam(parameters)
!
use contact_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "jeveux.h"
!
type(ContactParameters), intent(out) :: parameters
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get parameters from ContactComputaction.cxx/getData()
!
! --------------------------------------------------------------------------------------------------
!
    integer :: jcont
!
    call jevech('PCONFR', 'L', jcont)
!
! - Contact
!
    parameters%algo_cont = nint(zr(jcont+23))
    parameters%type_cont = nint(zr(jcont+24))
    parameters%vari_cont = zr(jcont+25)
    parameters%coef_cont = zr(jcont+26)
!
! - Friction
!
    parameters%l_fric = (zr(jcont+30) > 0.5d0)
    parameters%algo_fric = nint(zr(jcont+31))
    parameters%type_fric = nint(zr(jcont+32))
    parameters%coef_fric = zr(jcont+33)
    parameters%threshold_given = zr(jcont+34)
!
! - Other
!
    parameters%proj_tole = zr(jcont+40)
!
end subroutine
