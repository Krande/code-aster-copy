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
    function mefin5(nbz, nbgrp, imod, icyl, jmod,&
                    jcyl, z, f1, f2, f3,&
                    g)
        integer(kind=8) :: nbgrp
        integer(kind=8) :: nbz
        integer(kind=8) :: imod
        integer(kind=8) :: icyl
        integer(kind=8) :: jmod
        integer(kind=8) :: jcyl
        real(kind=8) :: z(*)
        real(kind=8) :: f1(nbz*nbgrp, *)
        real(kind=8) :: f2(nbz*nbgrp, *)
        real(kind=8) :: f3(*)
        real(kind=8) :: g(*)
        real(kind=8) :: mefin5
    end function mefin5
end interface
