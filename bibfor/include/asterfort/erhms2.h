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
    subroutine erhms2(ino, nbs, theta, jac,&
                      nx, ny, sielnp, adsip, sielnm,&
                      nbcmp, typmav, tbref1, tbref2, ivois,&
                      tm2h1s)
        integer(kind=8) :: ino
        integer(kind=8) :: nbs
        real(kind=8) :: theta
        real(kind=8) :: jac(3)
        real(kind=8) :: nx(3)
        real(kind=8) :: ny(3)
        real(kind=8) :: sielnp(140)
        integer(kind=8) :: adsip
        real(kind=8) :: sielnm(140)
        integer(kind=8) :: nbcmp
        character(len=8) :: typmav
        integer(kind=8) :: tbref1(12)
        integer(kind=8) :: tbref2(12)
        integer(kind=8) :: ivois
        real(kind=8) :: tm2h1s(3)
    end subroutine erhms2
end interface
