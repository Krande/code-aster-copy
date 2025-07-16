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
    subroutine pteequ(nume_equa    , base, neq, igds, nb_cmp_field,&
                      field_to_cata)
        character(len=19), intent(in) :: nume_equa
        integer(kind=8), intent(in) :: neq
        integer(kind=8), intent(in) :: igds
        integer(kind=8), intent(in) :: nb_cmp_field
        integer(kind=8), pointer :: field_to_cata(:)
        character(len=1), intent(in) :: base
    end subroutine pteequ
end interface
