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
    subroutine mefasr(ndim, nbcyl, nbgrp, nbtron, numgrp,&
                      idir, igrp, xint, yint, rint,&
                      sgn, orig, beta, a, b)
        integer(kind=8) :: nbtron
        integer(kind=8) :: nbcyl
        integer(kind=8) :: ndim(14)
        integer(kind=8) :: nbgrp
        integer(kind=8) :: numgrp(*)
        integer(kind=8) :: idir
        integer(kind=8) :: igrp
        real(kind=8) :: xint(*)
        real(kind=8) :: yint(*)
        real(kind=8) :: rint(*)
        integer(kind=8) :: sgn(*)
        integer(kind=8) :: orig(*)
        real(kind=8) :: beta(*)
        real(kind=8) :: a(2*nbtron*nbcyl, *)
        real(kind=8) :: b(*)
    end subroutine mefasr
end interface
