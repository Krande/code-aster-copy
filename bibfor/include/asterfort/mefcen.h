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
interface
    subroutine mefcen(caelem, iequiv, nbcyl, nbz, irot,&
                      numnog, nbnog, nummag, numgrp, coor,&
                      cent, req, xint, yint, zint,&
                      rint, nbgrp)
        integer(kind=8) :: nbgrp
        integer(kind=8) :: nbz
        integer(kind=8) :: nbcyl
        character(len=19) :: caelem
        integer(kind=8) :: iequiv
        integer(kind=8) :: irot(3)
        integer(kind=8) :: numnog(*)
        integer(kind=8) :: nbnog(*)
        integer(kind=8) :: nummag(*)
        integer(kind=8) :: numgrp(*)
        real(kind=8) :: coor(*)
        real(kind=8) :: cent(2*nbcyl)
        real(kind=8) :: req(nbgrp)
        real(kind=8) :: xint(nbcyl)
        real(kind=8) :: yint(nbcyl)
        real(kind=8) :: zint(nbz, nbgrp)
        real(kind=8) :: rint(nbcyl)
    end subroutine mefcen
end interface
