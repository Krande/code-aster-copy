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
    subroutine ccchel(option, &
                      modelZ, materFieldZ, materCodeZ, caraElemZ, listLoadZ, &
                      resultIn, resultOut, resultType, &
                      numeStore, numeStorePrev, &
                      ligrel, isTransient, postCompPoux, jvBase, &
                      fieldNameOut)
        use postComp_type
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ, listLoadZ
        character(len=8), intent(in) :: resultIn, resultOut
        character(len=16), intent(in) ::resultType
        integer(kind=8), intent(in) :: numeStore, numeStorePrev
        character(len=24), intent(in) :: ligrel
        aster_logical, intent(in) :: isTransient
        type(POST_COMP_POUX), intent(in) :: postCompPoux
        character(len=1), intent(in) :: jvBase
        character(len=24), intent(out) :: fieldNameOut
    end subroutine ccchel
end interface
