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
    subroutine ccchci(questz, type_comp, crit, norm, nb_form, &
                      repi)
        character(len=*), intent(in) :: questz
        character(len=16), intent(in) :: type_comp
        character(len=16), intent(in) :: crit
        character(len=16), intent(in) :: norm
        integer(kind=8), intent(in) :: nb_form
        integer(kind=8), intent(out) :: repi
    end subroutine ccchci
end interface
