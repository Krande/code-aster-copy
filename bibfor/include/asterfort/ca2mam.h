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
    subroutine ca2mam(modelInterfaceZ, incr, lchin, lpain, &
                      numeDof, matrAsse)
        character(len=*), intent(in) :: modelInterfaceZ
        character(len=3), intent(in) :: incr
        character(len=8), intent(in) :: lpain(2)
        character(len=24), intent(in) :: lchin(2)
        character(len=14), intent(out) :: numeDof
        character(len=24), intent(out) :: matrAsse
    end subroutine ca2mam
end interface
