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
    subroutine utflm2(mailla, tabmai, nbma, dim, typmai,&
                      nbtrou, tatrou)
        character(len=8), intent(in) :: mailla
        integer(kind=8), intent(in) :: nbma
        integer(kind=8), intent(in) :: tabmai(nbma)
        integer(kind=8), intent(in) :: dim
        character(len=*), intent(in) :: typmai
        integer(kind=8), intent(out) :: nbtrou
        integer(kind=8), intent(out) :: tatrou(nbma)
    end subroutine utflm2
end interface
