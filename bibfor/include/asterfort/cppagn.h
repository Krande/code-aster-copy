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
!
interface
    subroutine cppagn(main, maout, nbma, lima, izone, typ_dec, jcninv,&
                      same_zone, dec)
        character(len=8), intent(in) :: main
        character(len=8), intent(in) :: maout
        integer(kind=8), intent(in) :: nbma
        integer(kind=8), intent(in) :: lima(nbma)
        integer(kind=8), intent(in) :: izone
        integer(kind=8), intent(in) :: typ_dec
        integer(kind=8), intent(in) :: jcninv
        aster_logical, intent(in) :: same_zone
        integer(kind=8), intent(inout) :: dec
    end subroutine cppagn
end interface
