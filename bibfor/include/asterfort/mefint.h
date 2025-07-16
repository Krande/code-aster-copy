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
    subroutine mefint(nbz, nbgrp, nbmod, nbnoe, nbddl,&
                      irot, numnog, nbnog, zint, defm,&
                      phix, phiy, z, num)
        integer(kind=8) :: nbnoe
        integer(kind=8) :: nbmod
        integer(kind=8) :: nbgrp
        integer(kind=8) :: nbz
        integer(kind=8) :: nbddl
        integer(kind=8) :: irot(3)
        integer(kind=8) :: numnog(nbgrp)
        integer(kind=8) :: nbnog(nbgrp)
        real(kind=8) :: zint(nbz, nbgrp)
        real(kind=8) :: defm(6*nbnoe, nbmod)
        real(kind=8) :: phix(nbz, nbgrp, nbmod)
        real(kind=8) :: phiy(nbz, nbgrp, nbmod)
        real(kind=8) :: z(*)
        integer(kind=8) :: num(nbz)
    end subroutine mefint
end interface
