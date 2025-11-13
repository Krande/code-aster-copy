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
interface
    subroutine vdefge(nomte, nb1, npgsr, xr, epais, &
                      sigmElno, efgeElno)
        character(len=16), intent(in) :: nomte
        integer(kind=8), intent(in) :: nb1, npgsr
        real(kind=8), intent(in) :: xr(*), epais
        real(kind=8), intent(in) :: sigmElno(6, 27)
        real(kind=8), intent(out) :: efgeElno(8, 9)
    end subroutine vdefge
end interface
