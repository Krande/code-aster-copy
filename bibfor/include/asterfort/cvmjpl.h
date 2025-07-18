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
    subroutine cvmjpl(mod, nmat, mater, timed, timef,&
                      epsd, deps, sigf, vinf, sigd,&
                      vind, nvi, nr, dsde)
        common/tdim/ ndt,ndi
        integer(kind=8) :: ndt
        integer(kind=8) :: ndi
        integer(kind=8) :: nr
        integer(kind=8) :: nvi
        integer(kind=8) :: nmat
        character(len=8) :: mod
        real(kind=8) :: mater(nmat, 2)
        real(kind=8) :: timed
        real(kind=8) :: timef
        real(kind=8) :: epsd(*)
        real(kind=8) :: deps(*)
        real(kind=8) :: sigf(*)
        real(kind=8) :: vinf(*)
        real(kind=8) :: sigd(*)
        real(kind=8) :: vind(*)
        real(kind=8) :: dsde(6, 6)
    end subroutine cvmjpl
end interface
