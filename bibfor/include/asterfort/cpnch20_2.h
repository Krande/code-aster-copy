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
    subroutine cpnch20_2(main,numa,coor,ind,nomnoe, conneo)
        integer(kind=8), intent(in) :: ind
        integer(kind=8), intent(in) :: numa
        real(kind=8),intent(out) :: coor(3, *)
        character(len=8), intent(in) :: main
        character(len=24), intent(in) :: nomnoe
        integer(kind=8), intent(in) :: conneo(*)
    end subroutine cpnch20_2
end interface
