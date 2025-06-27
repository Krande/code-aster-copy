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
    subroutine ermes3(noe, ifa, tymvol, nnof, typmav,&
                      iref1, ivois, isig, nbcmp, dsg11,&
                      dsg22, dsg33, dsg12, dsg13, dsg23)
        integer(kind=8) :: noe(9, 6, 4)
        integer(kind=8) :: ifa
        integer(kind=8) :: tymvol
        integer(kind=8) :: nnof
        character(len=8) :: typmav
        integer(kind=8) :: iref1
        integer(kind=8) :: ivois
        integer(kind=8) :: isig
        integer(kind=8) :: nbcmp
        real(kind=8) :: dsg11(9)
        real(kind=8) :: dsg22(9)
        real(kind=8) :: dsg33(9)
        real(kind=8) :: dsg12(9)
        real(kind=8) :: dsg13(9)
        real(kind=8) :: dsg23(9)
    end subroutine ermes3
end interface
