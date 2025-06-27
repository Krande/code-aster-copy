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
    subroutine ermeb3(noe, ifa, tymvol, nnof, iref1,&
                      iref2, ivois, igeom, isig, nbcmp,&
                      inst, nx, ny, nz, sig11,&
                      sig22, sig33, sig12, sig13, sig23,&
                      chx, chy, chz)
        integer(kind=8) :: noe(9, 6, 4)
        integer(kind=8) :: ifa
        integer(kind=8) :: tymvol
        integer(kind=8) :: nnof
        integer(kind=8) :: iref1
        integer(kind=8) :: iref2
        integer(kind=8) :: ivois
        integer(kind=8) :: igeom
        integer(kind=8) :: isig
        integer(kind=8) :: nbcmp
        real(kind=8) :: inst
        real(kind=8) :: nx(9)
        real(kind=8) :: ny(9)
        real(kind=8) :: nz(9)
        real(kind=8) :: sig11(9)
        real(kind=8) :: sig22(9)
        real(kind=8) :: sig33(9)
        real(kind=8) :: sig12(9)
        real(kind=8) :: sig13(9)
        real(kind=8) :: sig23(9)
        real(kind=8) :: chx(9)
        real(kind=8) :: chy(9)
        real(kind=8) :: chz(9)
    end subroutine ermeb3
end interface
