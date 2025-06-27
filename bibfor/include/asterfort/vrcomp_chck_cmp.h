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
!
interface
    subroutine vrcomp_chck_cmp(mesh, nbCell,&
                               comporCurrZ,&
                               comporCurr, comporPrev,&
                               variRedu, comp_comb_2,&
                               ligrelCurr, ligrelPrev,&
                               verbose,&
                               nbSpgDifferent, nbVariDifferent, l_modif_vari)
        character(len=8), intent(in) :: mesh
        integer(kind=8), intent(in) :: nbCell
        character(len=*), intent(in) :: comporCurrZ
        character(len=19), intent(in) :: comporCurr, comporPrev
        character(len=19), intent(in) :: variRedu
        character(len=48), intent(in) :: comp_comb_2
        character(len=19), intent(in) :: ligrelCurr, ligrelPrev
        aster_logical, intent(in) :: verbose
        aster_logical, intent(out) :: nbSpgDifferent, nbVariDifferent, l_modif_vari
    end subroutine vrcomp_chck_cmp
end interface
