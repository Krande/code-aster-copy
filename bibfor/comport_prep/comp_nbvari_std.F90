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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_nbvari_std(rela_comp, defo_comp, type_cpla,&
                           kit_comp , post_iter, mult_comp,&
                           regu_visc,&
                           l_cristal, l_implex ,&
                           nb_vari  , nume_comp)
!
implicit none
!
#include "asterf_types.h"
#include "asterc/lcinfo.h"
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_code.h"
#include "asterfort/jeveuo.h"
!
character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla
character(len=16), intent(in) :: kit_comp(4), post_iter
character(len=16), intent(in) :: mult_comp, regu_visc
aster_logical, intent(in) :: l_cristal, l_implex
integer, intent(inout) :: nume_comp(4)
integer, intent(out) :: nb_vari
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get number of internal variables for standard constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  type_cpla        : plane stress method
! In  kit_comp         : KIT comportment
! In  post_iter        : type of post_treatment
! In  mult_comp        : multi-comportment
! In  l_cristal        : .true. if *CRISTAL comportment
! In  l_implex         : .true. if IMPLEX method
! IO  nume_comp        : number LCxxxx subroutine
! Out nb_vari          : number of internal variables
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: comp_code_py
    integer :: idummy
    character(len=8) :: sdcomp
    integer :: nb_vari_cris
    integer, pointer :: v_cpri(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_vari = 0

! - Coding composite comportment (Python)
    call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp,&
                        post_iter, regu_visc, l_implex,&
                        comp_code_py)

! - Get number of total internal state variables and index of law
    call lcinfo(comp_code_py, nume_comp(1), nb_vari, idummy)

! - Special for CRISTAL
    if (l_cristal) then
        sdcomp = mult_comp(1:8)
        call jeveuo(sdcomp//'.CPRI', 'L', vi=v_cpri)
        nb_vari_cris = v_cpri(3)
        nb_vari      = nb_vari + nb_vari_cris
        if (defo_comp .eq. 'SIMO_MIEHE') then
            nb_vari = nb_vari + 3 + 9
        endif
    endif

! - End of encoding
    call lcdiscard(comp_code_py)
!
end subroutine
