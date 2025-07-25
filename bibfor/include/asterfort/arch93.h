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
    subroutine arch93(resultName, concep, nume, raide, nbmodd,&
                      nbmodf, nbmoda, nbmoad, nbmodi, nbpsmo)
        character(len=8), intent(in) :: resultName
        character(len=16) :: concep
        character(len=14) :: nume
        character(len=19) :: raide
        integer(kind=8) :: nbmodd
        integer(kind=8) :: nbmodf
        integer(kind=8) :: nbmoda
        integer(kind=8) :: nbmoad
        integer(kind=8) :: nbmodi
        integer(kind=8) :: nbpsmo
    end subroutine arch93
end interface
