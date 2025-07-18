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
subroutine thm_kit_nvar(rela_thmc, rela_hydr, rela_meca, rela_ther, &
                        nb_vari_thmc, nb_vari_hydr, nb_vari_meca, nb_vari_ther, &
                        nume_comp_thmc, nume_comp_hydr, nume_comp_meca, nume_comp_ther)
!
    implicit none
!
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
!
    character(len=16), intent(in) :: rela_thmc
    character(len=16), intent(in) :: rela_hydr
    character(len=16), intent(in) :: rela_meca
    character(len=16), intent(in) :: rela_ther
    integer(kind=8), intent(out) :: nb_vari_thmc
    integer(kind=8), intent(out) :: nb_vari_hydr
    integer(kind=8), intent(out) :: nb_vari_meca
    integer(kind=8), intent(out) :: nb_vari_ther
    integer(kind=8), intent(out) :: nume_comp_thmc
    integer(kind=8), intent(out) :: nume_comp_hydr
    integer(kind=8), intent(out) :: nume_comp_meca
    integer(kind=8), intent(out) :: nume_comp_ther
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_thmc        : relation for coupling
! In  rela_hydr        : relation for hydraulic
! In  rela_meca        : relation for mechanic
! In  rela_ther        : relation for thermic
! Out nb_vari_thmc     : number of internal variables for coupling
! Out nb_vari_hydr     : number of internal variables for hydraulic
! Out nb_vari_meca     : number of internal variables for mechanic
! Out nb_vari_ther     : number of internal variables for thermic
! Out nume_comp_thmc   : number LCxxxx subroutine for coupling
! Out nume_comp_hydr   : number LCxxxx subroutine for hydraulic
! Out nume_comp_meca   : number LCxxxx subroutine for mechanic
! Out nume_comp_ther   : number LCxxxx subroutine for thermic
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: rela_thmc_py, rela_ther_py, rela_hydr_py, rela_meca_py
    integer(kind=8) :: ibid
!
! --------------------------------------------------------------------------------------------------
!
    nb_vari_thmc = 0
    nb_vari_ther = 0
    nb_vari_hydr = 0
    nb_vari_meca = 0
    nume_comp_thmc = 0
    nume_comp_ther = 0
    nume_comp_hydr = 0
    nume_comp_meca = 0
    call lccree(1, rela_thmc, rela_thmc_py)
    call lcinfo(rela_thmc_py, nume_comp_thmc, nb_vari_thmc, ibid)
    call lcdiscard(rela_thmc_py)
    call lccree(1, rela_ther, rela_ther_py)
    call lcinfo(rela_ther_py, nume_comp_ther, nb_vari_ther, ibid)
    call lcdiscard(rela_ther_py)
    call lccree(1, rela_hydr, rela_hydr_py)
    call lcinfo(rela_hydr_py, nume_comp_hydr, nb_vari_hydr, ibid)
    call lcdiscard(rela_hydr_py)
    call lccree(1, rela_meca, rela_meca_py)
    call lcinfo(rela_meca_py, nume_comp_meca, nb_vari_meca, ibid)
    call lcdiscard(rela_meca_py)
!
end subroutine
