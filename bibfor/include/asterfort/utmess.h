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
    subroutine utmess(typ, idmess, nk, valk, sk, &
                      ni, vali, si, nr, valr, &
                      sr, num_except, fname)
        character(len=*), intent(in) :: typ
        character(len=*), intent(in) :: idmess
        integer(kind=8), intent(in), optional :: nk
        character(len=*), intent(in), optional, target :: valk(*)
        character(len=*), intent(in), optional :: sk
        integer(kind=8), intent(in), optional :: ni
        integer(kind=8), intent(in), optional, target :: vali(*)
        integer(kind=8), intent(in), optional :: si
        integer(kind=8), intent(in), optional :: nr
        real(kind=8), intent(in), optional, target :: valr(*)
        real(kind=8), intent(in), optional :: sr
        integer(kind=8), intent(in), optional :: num_except
        character(len=*), optional :: fname
    end subroutine utmess
end interface
