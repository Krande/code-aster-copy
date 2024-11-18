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
    subroutine vefpme(stop, &
                      modelZ, caraElemZ, materFieldZ, matecoZ, &
                      loadNameJvZ, loadInfoJvZ, &
                      timePara, &
                      dispPrevZ, dispCumuInstZ, &
                      varcCurrZ, &
                      vectElemZ, ligrelCalcZ_, jvBase_)
        character(len=1), intent(in) :: stop
        character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, matecoZ
        character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
        real(kind=8), intent(in) :: timePara(3)
        character(len=*), intent(in) :: dispPrevZ, dispCumuInstZ
        character(len=*), intent(in) :: varcCurrZ
        character(len=*), intent(inout) :: vectElemZ
        character(len=*), optional, intent(in) :: ligrelCalcZ_
        character(len=1), optional, intent(in) :: jvBase_
    end subroutine vefpme
end interface
