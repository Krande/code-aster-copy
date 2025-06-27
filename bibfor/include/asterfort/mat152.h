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
    subroutine mat152(option, modelDime, modelInterface, ivalk, &
                      nbMode, matrAsseX, matrAsseY, matrAsseZ, numeDof)
        character(len=9), intent(in) :: option
        character(len=2), intent(in) :: modelDime
        character(len=8), intent(in) :: modelInterface
        integer(kind=8), intent(in) :: ivalk, nbMode
        character(len=14), intent(out) :: numeDof
        character(len=19), intent(out) :: matrAsseX, matrAsseY, matrAsseZ
    end subroutine mat152
end interface
