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
!
#include "asterf_types.h"
!
interface
    subroutine utin3d(igeom, nsomm, ino, ityp, inst,&
                      insold, k8cart, ltheta, niv, ifm,&
                      option, valfp, valfm, noe)
        integer(kind=8) :: igeom
        integer(kind=8) :: nsomm
        integer(kind=8) :: ino
        integer(kind=8) :: ityp
        real(kind=8) :: inst
        real(kind=8) :: insold
        character(len=8) :: k8cart
        aster_logical :: ltheta
        integer(kind=8) :: niv
        integer(kind=8) :: ifm
        integer(kind=8) :: option
        real(kind=8) :: valfp(9)
        real(kind=8) :: valfm(9)
        integer(kind=8) :: noe(9, 6, 3)
    end subroutine utin3d
end interface
