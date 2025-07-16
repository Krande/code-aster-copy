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
subroutine comp_nbvari_kit(kit_comp, &
                           l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg, &
                           nbVariFromKit, nbVariKit, numeLawKit)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cg_kit_nvar.h"
#include "asterfort/ddi_kit_nvar.h"
#include "asterfort/thm_kit_nvar.h"
#include "asterfort/meta_kit_nvar.h"
!
    character(len=16), intent(in) :: kit_comp(4)
    aster_logical, intent(in) :: l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg
    integer(kind=8), intent(out) :: nbVariFromKit, numeLawKit(4), nbVariKit(4)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Count number of internal variables for KIT
!
! --------------------------------------------------------------------------------------------------
!
! In  kit_comp         : KIT comportment
! In  l_kit_meta       : .true. if kit metallurgy
! In  l_kit_thm        : .true. if kit THM
! In  l_kit_ddi        : .true. if kit DDI
! In  l_kit_cg         : .true. if kit CG
! Out nbVariFromKit    : total number of internal state variables from kit
! Out numeLawKit       : index of subroutine for components in kit
! Out nbVariKit        : number of internal state variables for components in kit
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: rela_thmc, rela_ther, rela_hydr, rela_meca
    character(len=16) :: metaRela, metaPhas, metaGlob
    integer(kind=8) :: nbVariMetaRela, nbMetaPhas, nbVariMetaGlob
    integer(kind=8) :: nb_vari_thmc, nb_vari_ther, nb_vari_hydr, nb_vari_meca
    character(len=16) :: rela_flua, rela_plas, rela_cpla, rela_coup
    integer(kind=8) :: nb_vari_flua, nb_vari_plas, nb_vari_cpla, nb_vari_coup
    character(len=16) :: rela_comp_cg(2)
    integer(kind=8) :: nb_vari_cg(2), numeCompCG(2)
    integer(kind=8) :: nume_comp_plas, nume_comp_flua
    integer(kind=8) :: nume_comp_thmc, nume_comp_hydr, nume_comp_meca, nume_comp_ther
!
! --------------------------------------------------------------------------------------------------
!
    nbVariFromKit = 0
    nbVariKit = 0
    numeLawKit = 0

! - Number of internal state variables for KIT THM
    if (l_kit_thm) then
        rela_meca = kit_comp(1)
        rela_hydr = kit_comp(2)
        rela_ther = kit_comp(3)
        rela_thmc = kit_comp(4)
        call thm_kit_nvar(rela_thmc, rela_hydr, rela_meca, rela_ther, &
                          nb_vari_thmc, nb_vari_hydr, nb_vari_meca, nb_vari_ther, &
                          nume_comp_thmc, nume_comp_hydr, nume_comp_meca, nume_comp_ther)
        nbVariKit(1) = nb_vari_thmc
        nbVariKit(2) = nb_vari_ther
        nbVariKit(3) = nb_vari_hydr
        nbVariKit(4) = nb_vari_meca
        numeLawKit(1) = nume_comp_thmc
        numeLawKit(2) = nume_comp_ther
        numeLawKit(3) = nume_comp_hydr
        numeLawKit(4) = nume_comp_meca
    end if

! - Number of internal state variables for KIT META
    if (l_kit_meta) then
        metaPhas = kit_comp(1)
        metaRela = kit_comp(2)
        metaGlob = kit_comp(3)
        call meta_kit_nvar(metaPhas, metaRela, metaGlob, &
                           nbMetaPhas, nbVariMetaRela, nbVariMetaGlob)
        nbVariKit(1) = nbVariMetaRela
        nbVariKit(2) = nbMetaPhas
        nbVariKit(3) = nbVariMetaGlob
        nbVariFromKit = nbVariMetaGlob+nbVariMetaRela*nbMetaPhas
    end if

! - Number of internal state variables for KIT DDI
    if (l_kit_ddi) then
        rela_flua = kit_comp(1)
        rela_plas = kit_comp(2)
        rela_cpla = kit_comp(3)
        rela_coup = kit_comp(4)
        call ddi_kit_nvar(rela_flua, rela_plas, rela_cpla, rela_coup, &
                          nb_vari_flua, nb_vari_plas, nb_vari_cpla, nb_vari_coup, &
                          nume_comp_plas, nume_comp_flua)
        nbVariKit(1) = nb_vari_flua
        nbVariKit(2) = nb_vari_plas
        nbVariKit(3) = nb_vari_cpla
        nbVariKit(4) = nb_vari_coup
        numeLawKit(1) = nume_comp_flua
        numeLawKit(2) = nume_comp_plas
    end if

! - Number of internal state variables for KIT CG
    if (l_kit_cg) then
        rela_comp_cg(1) = kit_comp(1)
        rela_comp_cg(2) = kit_comp(2)
        call cg_kit_nvar(rela_comp_cg, nb_vari_cg, numeCompCG)
        nbVariKit(1) = nb_vari_cg(1)
        nbVariKit(2) = nb_vari_cg(2)
        numeLawKit(1) = numeCompCG(1)
        numeLawKit(2) = numeCompCG(2)
    end if
!
end subroutine
