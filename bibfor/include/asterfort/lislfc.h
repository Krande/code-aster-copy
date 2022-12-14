! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
    subroutine lislfc(list_load_resu, i_load      , i_excit   , l_load_user,&
                      l_func_c      , load_keyword, const_func, load_func, basez)
        character(len=19), intent(in) :: list_load_resu
        integer, intent(in) :: i_load, i_excit
        aster_logical, intent(in) :: l_load_user, l_func_c
        character(len=16), intent(in) :: load_keyword
        character(len=8), intent(inout) :: const_func
        character(len=8), intent(out) :: load_func
        character(len=1), intent(in), optional :: basez
    end subroutine lislfc
end interface
