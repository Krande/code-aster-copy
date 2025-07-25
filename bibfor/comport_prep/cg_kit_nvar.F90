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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine cg_kit_nvar(rela_comp_cg, nb_vari_cg, numeCompCG)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
!
    character(len=16), intent(in) :: rela_comp_cg(2)
    integer(kind=8), intent(out) :: nb_vari_cg(2)
    integer(kind=8), intent(out) :: numeCompCG(2)
!
! --------------------------------------------------------------------------------------------------
!
! KIT_CG
!
! Number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp_cg     : relations for KIT_CG
! Out nb_vari_cg       : number of internal variables for KIT_CG
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: rela_py
    integer(kind=8) :: numlc, nbvari_ext
!
! --------------------------------------------------------------------------------------------------
!
    nb_vari_cg = 0
    numeCompCG = 0
    call lccree(1, rela_comp_cg(1), rela_py)
    call lcinfo(rela_py, numeCompCG(1), nb_vari_cg(1), nbvari_ext)
    call lcdiscard(rela_py)
    call lccree(1, rela_comp_cg(2), rela_py)
    call lcinfo(rela_py, numeCompCG(2), nb_vari_cg(2), nbvari_ext)
    call lcdiscard(rela_py)
!
end subroutine
