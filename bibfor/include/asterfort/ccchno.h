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
#include "asterf_types.h"
!
interface
    subroutine ccchno(option, numeStore, resultIn, resultOut, fieldNameOut, &
                      lRestCell, nbRestCell, restCellJv, &
                      mesh, model, caraElem, optionBase, &
                      ligrel, ligrelHasBeenChanged, codret, &
                      fieldNameIn, ideb_, ifin_, vcham_)
        character(len=16), intent(in) :: option
        integer(kind=8), intent(in) :: numeStore
        character(len=8), intent(in) :: resultIn, resultOut
        character(len=24), intent(out) :: fieldNameOut
        aster_logical, intent(in) :: lRestCell
        integer(kind=8), intent(in) :: nbRestCell
        character(len=24), intent(in) :: restCellJv
        character(len=8), intent(in) :: mesh, model, caraElem
        character(len=1), intent(in) :: optionBase
        character(len=24), intent(in) :: ligrel
        integer(kind=8), intent(out) :: codret
        aster_logical, intent(in) :: ligrelHasBeenChanged
        character(len=19), intent(in) :: fieldNameIn
        integer(kind=8), optional, intent(in) :: ideb_, ifin_
        character(len=24), optional, intent(in) :: vcham_
    end subroutine ccchno
end interface
