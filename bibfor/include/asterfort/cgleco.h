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
!
#include "asterf_types.h"
!
interface
    subroutine cgleco(resu, modele, mate, iord0, compor,incr)
        integer(kind=8), intent(in) :: iord0
        character(len=8), intent(in) :: resu
        character(len=8), intent(in) :: modele
        character(len=8), intent(in) :: mate
        character(len=19), intent(out) :: compor
        aster_logical, intent(out) :: incr
    end subroutine cgleco
end interface
