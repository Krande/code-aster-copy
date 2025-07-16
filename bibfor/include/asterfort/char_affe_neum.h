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
    subroutine char_affe_neum(model      , mesh, ndim ,&
                              keywordfact, iocc,&
                              nbMap      , map , nbCmp)
        character(len=8), intent(in) :: model, mesh
        integer(kind=8), intent(in) :: ndim
        character(len=16), intent(in) :: keywordfact
        integer(kind=8), intent(in) :: iocc, nbMap
        character(len=19), intent(in) :: map(nbMap)
        integer(kind=8), intent(in) :: nbCmp(nbMap)
    end subroutine char_affe_neum
end interface
