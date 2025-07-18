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
subroutine ddi_kit_nvar(rela_flua, rela_plas, rela_cpla, rela_coup, &
                        nb_vari_flua, nb_vari_plas, nb_vari_cpla, nb_vari_coup, &
                        nume_comp_plas, nume_comp_flua)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
!
    character(len=16), intent(in) :: rela_flua
    character(len=16), intent(in) :: rela_plas
    character(len=16), intent(in) :: rela_cpla
    character(len=16), intent(in) :: rela_coup
    integer(kind=8), intent(out) :: nb_vari_flua
    integer(kind=8), intent(out) :: nb_vari_plas
    integer(kind=8), intent(out) :: nb_vari_cpla
    integer(kind=8), intent(out) :: nb_vari_coup
    integer(kind=8), intent(out) :: nume_comp_plas
    integer(kind=8), intent(out) :: nume_comp_flua
!
! --------------------------------------------------------------------------------------------------
!
! KIT_DDI
!
! Number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_flua        : relation for creeping
! In  rela_plas        : relation for plasticity
! In  rela_cpla        : relation for plane stress (GLRC)
! In  rela_coup        : relation for coupling (GLRC)
! Out nb_vari_flua     : number of internal variables for creeping
! Out nb_vari_plas     : number of internal variables for plasticity
! Out nb_vari_cpla     : number of internal variables for plane stress
! Out nb_vari_coup     : number of internal variables for coupling
! Out nume_comp_plas   : number LCxxxx subroutine for plasticity
! Out nume_comp_flua   : number LCxxxx subroutine for creep
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: rela_py
    integer(kind=8) :: ibid, ibid2
!
! --------------------------------------------------------------------------------------------------
!
    nb_vari_flua = 0
    nb_vari_plas = 0
    nb_vari_cpla = 0
    nb_vari_coup = 0
    nume_comp_plas = 0
    nume_comp_flua = 0
    if (rela_flua .ne. ' ') then
        call lccree(1, rela_flua, rela_py)
        call lcinfo(rela_py, nume_comp_flua, nb_vari_flua, ibid)
        call lcdiscard(rela_py)
    end if
    if (rela_plas .ne. ' ') then
        call lccree(1, rela_plas, rela_py)
        call lcinfo(rela_py, nume_comp_plas, nb_vari_plas, ibid)
        call lcdiscard(rela_py)
    end if
    if (rela_cpla .ne. ' ') then
        call lccree(1, rela_cpla, rela_py)
        call lcinfo(rela_py, ibid, nb_vari_cpla, ibid2)
        call lcdiscard(rela_py)
    end if
    if (rela_coup .ne. ' ') then
        call lccree(1, rela_coup, rela_py)
        call lcinfo(rela_py, ibid, nb_vari_coup, ibid2)
        call lcdiscard(rela_py)
    end if
!
end subroutine
