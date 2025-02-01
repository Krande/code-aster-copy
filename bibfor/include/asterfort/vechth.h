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
    subroutine vechth(typeTher, &
                      modelZ, matecoZ, &
                      loadNameJvZ, loadInfoJvZ, &
                      timeCurr, &
                      vectElemZ, &
                      varcCurrZ_, timeMapZ_, tempPrevZ_, timeMoveZ_, &
                      jvBase_)
        character(len=4), intent(in) :: typeTher
        character(len=*), intent(in) :: modelZ, matecoZ
        character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
        real(kind=8), intent(in) :: timeCurr
        character(len=*), intent(inout) :: vectElemZ
        character(len=*), optional, intent(in) :: varcCurrZ_, timeMapZ_, tempPrevZ_, timeMoveZ_
        character(len=1), optional, intent(in) :: jvBase_
    end subroutine vechth
end interface
