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
subroutine comp_nbvari(rela_comp, defo_comp, type_cpla, kit_comp, &
                       post_iter, mult_comp, regu_visc, post_incr, &
                       extern_type, extern_addr, &
                       nbVariUMAT, &
                       nbVari, numeLaw, nbVariKit, numeLawKit)
!
    implicit none
!
#include "asterc/mgis_get_sizeof_isvs.h"
#include "asterf_types.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_nbvari_kit.h"
#include "asterfort/comp_nbvari_std.h"
#include "asterfort/jeveuo.h"
!
    character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla
    character(len=16), intent(in) :: kit_comp(4), post_iter
    character(len=16), intent(in) :: mult_comp, regu_visc, post_incr, extern_addr
    integer(kind=8), intent(in) :: extern_type
    integer(kind=8), intent(in) :: nbVariUMAT
    integer(kind=8), intent(out) :: nbVari, numeLaw, nbVariKit(4), numeLawKit(4)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count the number of internal state variables and index of behaviours
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  type_cpla        : plane stress method
! In  kit_comp         : KIT comportment
! In  mult_comp        : multi-comportment (for crystal)
! In  post_iter        : type of post-treatment at each Newton iteration
! In  regu_visc        : keyword for viscuous regularization
! In  post_incr        : type of post-treatment at end of time step
! In  external_type    : type of type of integration (internal, official, proto, umat)
! In  external_ptr     : address of external behaviour
! In  nbVariUMAT       : number of internal state variables for UMAT
! Out nbVari           : number of internal state variables
! Out numeLaw          : index of subroutine for behaviour
! Out nbVariKit        : number of internal state variables for components in kit
! Out numeLawKit       : index of subroutine for components in kit
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbVariExte, nbVariFromKit, nbVariCrystal
    aster_logical :: l_cristal, l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg, l_kit
    integer(kind=8), pointer :: cpri(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = 0
    numeLaw = 0
    nbVariKit = 0
    numeLawKit = 0

! - Detection of specific cases
    call comp_meca_l(rela_comp, 'KIT', l_kit)
    call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)
    call comp_meca_l(rela_comp, 'KIT_META', l_kit_meta)
    call comp_meca_l(rela_comp, 'KIT_THM', l_kit_thm)
    call comp_meca_l(rela_comp, 'KIT_DDI', l_kit_ddi)
    call comp_meca_l(rela_comp, 'KIT_CG', l_kit_cg)

! - Get number of internal state variables for KIT
    nbVariFromKit = 0
    if (l_kit) then
        call comp_nbvari_kit(kit_comp, &
                             l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg, &
                             nbVariFromKit, nbVariKit, numeLawKit)
    end if

! - Special for CRISTAL
    nbVariCrystal = 0
    if (l_cristal) then
        call jeveuo(mult_comp(1:8)//'.CPRI', 'L', vi=cpri)
        nbVariCrystal = cpri(3)
        if (defo_comp .eq. 'SIMO_MIEHE') then
            nbVariCrystal = nbVariCrystal+3+9
        end if
    end if

! - Get number of internal state variables
    call comp_nbvari_std(rela_comp, defo_comp, type_cpla, &
                         kit_comp, post_iter, regu_visc, post_incr, &
                         nbVari, numeLaw)

! - Get number of internal state variables for external behaviours
    nbVariExte = 0
    if (extern_type .ne. 0) then
        if (extern_type .eq. 1 .or. extern_type .eq. 2) then
            ! MFront offi/proto
            call mgis_get_sizeof_isvs(extern_addr, nbVariExte)
            if (nbVariExte .eq. 0) then
                nbVariExte = 1
            end if
        else
            ! UMAT
            nbVariExte = nbVariUMAT
        end if
        nbVariKit(4) = nbVariExte
    end if

! - Total number of internal state variables
    nbVari = nbVariFromKit+nbVari
    nbVari = nbVariCrystal+nbVari
    nbVari = nbVariExte+nbVari
!
end subroutine
