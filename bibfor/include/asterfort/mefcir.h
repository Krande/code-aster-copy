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
    subroutine mefcir(ndim, nbcyl, nbgrp, numgrp, som,&
                      rint, dcent, ficent, d, fi,&
                      ppxx, ppxy, ppyx, ppyy, vnxx,&
                      vnxy, vnyx, vnyy, tmp)
        integer(kind=8) :: nbgrp
        integer(kind=8) :: nbcyl
        integer(kind=8) :: ndim(14)
        integer(kind=8) :: numgrp(*)
        real(kind=8) :: som(9)
        real(kind=8) :: rint(*)
        real(kind=8) :: dcent(*)
        real(kind=8) :: ficent(*)
        real(kind=8) :: d(*)
        real(kind=8) :: fi(*)
        real(kind=8) :: ppxx(nbcyl, nbgrp)
        real(kind=8) :: ppxy(nbcyl, nbgrp)
        real(kind=8) :: ppyx(nbcyl, nbgrp)
        real(kind=8) :: ppyy(nbcyl, nbgrp)
        real(kind=8) :: vnxx(nbcyl, nbgrp)
        real(kind=8) :: vnxy(nbcyl, nbgrp)
        real(kind=8) :: vnyx(nbcyl, nbgrp)
        real(kind=8) :: vnyy(nbcyl, nbgrp)
        real(kind=8) :: tmp(4, *)
    end subroutine mefcir
end interface
