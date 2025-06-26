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
    subroutine mecgme(stop, &
                      modelZ, caraElemZ, materFieldZ, matecoZ, comporZ, &
                      listLoadZ, &
                      timePrev, timeCurr, &
                      dispPrevZ, dispCumuInstZ, &
                      matrElemZ, &
                      varcCurrZ_, nharm_, &
                      ligrelCalcZ_, jvBase_)
        character(len=1), intent(in) :: stop
        character(len=*), intent(in) :: modelZ, caraElemZ
        character(len=*), intent(in) :: materFieldZ, matecoZ, comporZ
        character(len=*), intent(in) :: listLoadZ
        real(kind=8), intent(in) :: timePrev, timeCurr
        character(len=*), intent(in) :: dispPrevZ, dispCumuInstZ
        character(len=*), intent(in) :: matrElemZ
        character(len=*), optional, intent(in) :: varcCurrZ_
        integer(kind=8), optional, intent(in) :: nharm_
        character(len=*), optional, intent(in) :: ligrelCalcZ_
        character(len=1), optional, intent(in) :: jvBase_
    end subroutine mecgme
end interface
