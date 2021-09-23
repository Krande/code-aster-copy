! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

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
    subroutine comp_nbvari_std(rela_comp, defo_comp, type_cpla,&
                               kit_comp , post_iter, mult_comp,&
                               regu_visc,&
                               l_cristal, l_implex ,&
                               nb_vari  , nume_comp)
        character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla
        character(len=16), intent(in) :: kit_comp(4), post_iter
        character(len=16), intent(in) :: mult_comp, regu_visc
        aster_logical, intent(in) :: l_cristal, l_implex
        integer, intent(inout) :: nume_comp(4)
        integer, intent(out) :: nb_vari
    end subroutine comp_nbvari_std
end interface
