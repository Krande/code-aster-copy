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
subroutine nocart_c(carte, code, ncmp, groupma, mode, nma,&
                    limano, limanu, ligrel)

    character(len=*), intent(in) :: carte
    integer(kind=8), intent(in) :: code
    integer(kind=8), intent(in) :: ncmp
    character(len=*), intent(in) :: groupma
    character(len=*),intent(in) :: mode
    integer(kind=8), intent(in) :: nma
    character(len=*), intent(in) :: limano(*)
    integer(kind=8), intent(in) :: limanu(*)
    character(len=*), intent(in) ::  ligrel

end subroutine nocart_c
end interface
