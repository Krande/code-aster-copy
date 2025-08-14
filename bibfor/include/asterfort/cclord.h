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
    subroutine cclord(numeOptEff, nbStore, listStore, jvBaseName, isOptionFromUser, &
                      numeStoreMin, numeStoreMax, resultIn, resultOut)
        integer(kind=8), intent(in) :: numeOptEff
        integer(kind=8), intent(in) :: nbStore
        integer(kind=8), pointer :: listStore(:)
        character(len=8), intent(in) :: jvBaseName
        aster_logical, intent(in) :: isOptionFromUser
        integer(kind=8), intent(in) ::  numeStoreMin, numeStoreMax
        character(len=8), intent(in) :: resultIn, resultOut
    end subroutine cclord
end interface
