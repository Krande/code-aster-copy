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
#include "asterf_types.h"
interface
    subroutine merxth(l_stat, &
                      modelZ, caraElemZ, matecoZ, &
                      loadNameJvZ, loadInfoJvZ, &
                      tpsthe, timeMapZ, &
                      tempIterZ, comporTherZ, varcCurrZ, dryCurrZ, &
                      matrElemZ, jvBase)
        aster_logical, intent(in) :: l_stat
        character(len=*), intent(in) :: modelZ, caraElemZ, matecoZ
        character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
        real(kind=8), intent(in) :: tpsthe(6)
        character(len=*), intent(in) :: comporTherZ, timeMapZ
        character(len=*), intent(in) :: tempIterZ, varcCurrZ, dryCurrZ
        character(len=*), intent(in) :: matrElemZ
        character(len=1), intent(in) :: jvBase
    end subroutine merxth
end interface
