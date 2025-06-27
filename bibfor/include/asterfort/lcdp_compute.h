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

#include "asterf_types.h"
!
interface
    function lcdp_compute(resi,rigi,elas, itemax, prec, m, eps, ep, ka, state, &
                            s,deps_s, vip) &
        result(iret)

        use lcdp_module,         only: dp_material

        aster_logical,intent(in)     :: resi,rigi,elas
        type(dp_material),intent(in) :: m
        integer(kind=8),intent(in)           :: itemax
        real(kind=8),intent(in)      :: prec
        real(kind=8),intent(in)      :: eps(:)
        real(kind=8),intent(inout)   :: ep(:),ka
        integer(kind=8),intent(inout)        :: state
        real(kind=8),intent(out)     :: s(:),deps_s(:,:)
        integer(kind=8)                      :: iret
        real(kind=8) :: vip(9)
    end function
end interface
