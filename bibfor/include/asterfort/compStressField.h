! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
!
interface
    subroutine compStressField(result, model, mater, mateco, cara_elem, list_load, &
                           l_sief_elga, l_strx_elga, nbrank, times)
        character(len=*), intent(in) :: model, cara_elem, list_load, result, mater, mateco
        aster_logical, intent(in) :: l_sief_elga, l_strx_elga
        integer, intent(in) :: nbrank
        real(kind=8), intent(in) :: times(*)
    end subroutine compStressField
end interface
