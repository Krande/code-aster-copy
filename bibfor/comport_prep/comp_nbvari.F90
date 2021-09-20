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
subroutine comp_nbvari(rela_comp, defo_comp, type_cpla, kit_comp ,&
                       post_iter, meca_comp, mult_comp, regu_visc,&
                       l_implex ,&
                       libr_name, subr_name, model_dim, model_mfront,&
                       nbVariUMAT,&
                       nbVari, numeLaw, nbVariKit, numeLawKit)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_nbvari_std.h"
#include "asterfort/comp_nbvari_kit.h"
#include "asterfort/comp_nbvari_ext.h"
#include "asterfort/jeveuo.h"
!
character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla
character(len=16), intent(in) :: kit_comp(4), post_iter, meca_comp
character(len=16), intent(in) :: mult_comp, regu_visc
aster_logical, intent(in) :: l_implex
character(len=255), intent(in) :: libr_name, subr_name
integer, intent(in) :: model_dim
character(len=16), intent(in) :: model_mfront
integer, intent(in) :: nbVariUMAT
integer, intent(out) :: nbVari, numeLaw, nbVariKit(4), numeLawKit(4)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Count the number of internal state variables and index of behaviours
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  type_cpla        : plane stress method
! In  kit_comp         : KIT comportment
! In  post_iter        : type of post_treatment
! In  regu_visc        : keyword for viscuous regularization
! In  l_implex         : .true. if IMPLEX method
! In  mult_comp        : multi-comportment (for crystal)
! In  nbVariUMAT       : number of internal state variables for UMAT
! In  libr_name        : name of library if UMAT or MFront
! In  subr_name        : name of comportement in library if UMAT or MFront
! In  model_dim        : dimension of modelisation (2D or 3D)
! In  model_mfront     : type of modelisation MFront
! In  l_implex         : .true. if IMPLEX method
! Out nbVari           : number of internal state variables
! Out numeLaw          : index of subroutine for behaviour
! Out nbVariKit        : number of internal state variables for components in kit
! Out numeLawKit       : index of subroutine for components in kit
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbVariExte, nbVariFromKit, nbVariCrystal
    aster_logical :: l_cristal, l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg, l_kit
    aster_logical :: l_exte_comp, l_mfront_proto, l_mfront_offi, l_umat
    integer, pointer :: cpri(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = 0
    numeLaw = 0
    nbVariKit = 0
    numeLawKit = 0

! - Detection of specific cases
    call comp_meca_l(rela_comp, 'KIT'     , l_kit)
    call comp_meca_l(rela_comp, 'CRISTAL' , l_cristal)
    call comp_meca_l(rela_comp, 'KIT_META', l_kit_meta)
    call comp_meca_l(rela_comp, 'KIT_THM' , l_kit_thm)
    call comp_meca_l(rela_comp, 'KIT_DDI' , l_kit_ddi)
    call comp_meca_l(rela_comp, 'KIT_CG'  , l_kit_cg)

! - Get number of internal state variables for KIT
    nbVariFromKit = 0
    if (l_kit) then
        call comp_nbvari_kit(kit_comp,&
                             l_kit_meta   , l_kit_thm   , l_kit_ddi, l_kit_cg,&
                             nbVariFromKit, nbVariKit, numeLawKit)
    endif

! - Special for CRISTAL
    nbVariCrystal = 0
    if (l_cristal) then
        call jeveuo(mult_comp(1:8)//'.CPRI', 'L', vi=cpri)
        nbVariCrystal = cpri(3)
        if (defo_comp .eq. 'SIMO_MIEHE') then
            nbVariCrystal = nbVariCrystal + 3 + 9
        endif
    endif

! - Get number of internal state variables
    call comp_nbvari_std(rela_comp, defo_comp, type_cpla,&
                         kit_comp , post_iter, regu_visc,&
                         l_implex , nbVari   , numeLaw)

! - Get number of internal state variables for external behaviours
    nbVariExte = 0
    call comp_meca_l(meca_comp, 'EXTE_COMP'   , l_exte_comp)
    call comp_meca_l(meca_comp, 'MFRONT_PROTO', l_mfront_proto)
    call comp_meca_l(meca_comp, 'MFRONT_OFFI' , l_mfront_offi)
    call comp_meca_l(meca_comp, 'UMAT'        , l_umat)
    if (l_exte_comp) then
        call comp_nbvari_ext(l_umat        , nbVariUMAT   ,&
                             l_mfront_proto, l_mfront_offi,&
                             libr_name     , subr_name    ,&
                             model_dim     , model_mfront ,&
                             nbVariExte)
        nbVariKit(4) = nbVariExte
    endif

! - Total number of internal state variables
    nbVari = nbVariFromKit + nbVari
    nbVari = nbVariCrystal + nbVari
    nbVari = nbVariExte + nbVari
!
end subroutine
