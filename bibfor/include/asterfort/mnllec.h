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
    subroutine mnllec(imat, numedd, ordman, epsman, pasman, epscor, h,&
                      hf, itemax, nbran, nextr, epsbif)
        integer(kind=8) :: imat(2)
        character(len=24) :: numedd
        integer(kind=8) :: ordman
        real(kind=8) :: epsman
        integer(kind=8) :: pasman
        real(kind=8) :: epscor
        integer(kind=8) :: h
        integer(kind=8) :: hf
        integer(kind=8) :: itemax
        integer(kind=8) :: nbran
        integer(kind=8) :: nextr
        real(kind=8) :: epsbif
    end subroutine mnllec
end interface 
