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
    subroutine hujiid(mod, mater, indi, deps, i1e,&
                      yd, vind, dy, loop, dsig,&
                      bnews, mtrac, iret)
        character(len=8) :: mod
        real(kind=8) :: mater(22, 2)
        integer(kind=8) :: indi(7)
        real(kind=8) :: deps(6)
        real(kind=8) :: i1e
        real(kind=8) :: yd(18)
        real(kind=8) :: vind(*)
        real(kind=8) :: dy(18)
        aster_logical :: loop
        real(kind=8) :: dsig(6)
        aster_logical :: bnews(3)
        aster_logical :: mtrac
        integer(kind=8) :: iret
    end subroutine hujiid
end interface
