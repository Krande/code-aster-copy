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
    subroutine dstsie(option, fami, xyzl, pgl, depl,&
                      nbcou, cdl)
        character(len=16) :: option
        character(len=4) :: fami
        real(kind=8) :: xyzl(3, *)
        real(kind=8) :: pgl(3, *)
        real(kind=8) :: depl(*)
        integer(kind=8) :: nbcou
        real(kind=8) :: cdl(*)
    end subroutine dstsie
end interface
