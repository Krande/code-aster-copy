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
    subroutine zzcalb(igr, iel, npg, nno, wi,&
                      desc, sig, x, y, xmin,&
                      xmax, ymin, ymax, f)
        integer(kind=8) :: igr
        integer(kind=8) :: iel
        integer(kind=8) :: npg
        integer(kind=8) :: nno
        real(kind=8) :: wi(1)
        integer(kind=8) :: desc(1)
        real(kind=8) :: sig(1)
        real(kind=8) :: x(1)
        real(kind=8) :: y(1)
        real(kind=8) :: xmin
        real(kind=8) :: xmax
        real(kind=8) :: ymin
        real(kind=8) :: ymax
        real(kind=8) :: f(9, 4)
    end subroutine zzcalb
end interface
