! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
    subroutine nmdini(keywf  , list_inst     , tole,&
                      nb_inst, l_init_noexist, nume_ini )
        character(len=16), intent(in) :: keywf
        character(len=19), intent(in) :: list_inst
        real(kind=8), intent(in) :: tole
        integer, intent(in) :: nb_inst
        aster_logical, intent(out) :: l_init_noexist
        integer, intent(out) :: nume_ini
    end subroutine nmdini
end interface
