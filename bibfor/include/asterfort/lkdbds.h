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
    subroutine lkdbds(nmat, mater, i1, devsig, nvi,&
                      vint, para, val, dbetds, dbetdi)
        integer(kind=8) :: nvi
        integer(kind=8) :: nmat
        real(kind=8) :: mater(nmat, 2)
        real(kind=8) :: i1
        real(kind=8) :: devsig(6)
        real(kind=8) :: vint(nvi)
        real(kind=8) :: para(3)
        integer(kind=8) :: val
        real(kind=8) :: dbetds(6)
        real(kind=8) :: dbetdi
    end subroutine lkdbds
end interface
