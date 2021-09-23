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
subroutine comp_meca_name(nbVari     , nbVariMeca,&
                          l_excl     , vari_excl,&
                          l_kit_meta , l_mfront_offi, l_prot_comp,&
                          rela_comp  , defo_comp    , kit_comp     ,&
                          type_cpla  , post_iter    , regu_visc    ,&
                          libr_name  , subr_name    , model_mfront , model_dim,&
                          infoVari)
!
implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcvari.h"
#include "asterc/lcdiscard.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/comp_mfront_vname.h"
#include "asterfort/comp_meca_code.h"
!
integer, intent(in) :: nbVari, nbVariMeca
aster_logical, intent(in) :: l_excl
character(len=16), intent(in) :: vari_excl
aster_logical, intent(in) :: l_kit_meta, l_mfront_offi, l_prot_comp
character(len=16), intent(in) :: rela_comp, defo_comp, kit_comp(4)
character(len=16), intent(in) :: type_cpla, post_iter, regu_visc
character(len=255), intent(in) :: libr_name, subr_name
character(len=16), intent(in) :: model_mfront
integer, intent(in) :: model_dim
character(len=16), pointer :: infoVari(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  nbVari           : number of internal variables
! In  nbVariMeca       : number of internal variables for mechanic
! In  l_excl           : .true. if exception case (no names for internal variables)
! In  vari_excl        : name of internal variables if l_excl
! In  l_kit_meta       : .true. if metallurgy
! In  l_mfront_offi    : .true. if MFront official
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  kit_comp         : KIT comportment
! In  type_cpla        : plane stress method
! In  post_iter        : type of post_treatment
! In  regu_visc        : keyword for viscuous regularization
! In  libr_name        : name of library
! In  subr_name        : name of comportement in library
! In  model_mfront     : type of modelisation MFront
! In  model_dim        : dimension of modelisation (2D or 3D)
! Ptr infoVari         : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! - No internal state variable for IMPLEX
    aster_logical, parameter :: l_implex = ASTER_FALSE
    character(len=6) :: metaPhasName(10)
    character(len=8) :: metaRelaName(30)
    character(len=16) :: metaGlobName(30)
    integer :: idummy, idummy2, nbVariOther, iVariMeca, iVari
    character(len=16) :: compCodePy
    character(len=16) :: metaPhas, metaRela, metaGlob
    character(len=16) :: metaPhasPy, metaRelaPy, metaGlobPy
    integer :: nbMetaPhas, nbVariMetaRela, nbVariMetaGlob
    integer :: iMetaPhas, iVariMetaRela, iVariMetaGlob
!
! --------------------------------------------------------------------------------------------------
!
    if (l_excl) then
        infoVari(1:nbVari) = vari_excl
    else
! ----- Name of internal state variables
        if (l_prot_comp) then
            call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp,&
                                post_iter, regu_visc, l_implex,&
                                compCodePy)
            nbVariOther = nbVari - nbVariMeca
            do iVariMeca = 1, nbVariMeca
                infoVari(iVariMeca) = 'NoName'
            end do
            if (nbVariOther .ne. 0) then
                call lcvari(compCodePy, nbVariOther, infoVari(nbVariMeca+1:nbVari))
            endif
            call lcdiscard(compCodePy)

        else if (l_mfront_offi) then
            call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp,&
                                post_iter, regu_visc, l_implex,&
                                compCodePy)
            nbVariOther = nbVari - nbVariMeca
            call comp_mfront_vname(nbVariMeca, &
                                   libr_name , subr_name, model_mfront, model_dim,&
                                   infoVari)
            if (nbVariOther .ne. 0) then
                call lcvari(compCodePy, nbVariOther, infoVari(nbVariMeca+1:nbVari))
            endif
            call lcdiscard(compCodePy)

        else
            if (l_kit_meta) then
                metaPhas = kit_comp(1)
                metaRela = kit_comp(2)
                metaGlob = kit_comp(3)
                call lccree(1, metaPhas, metaPhasPy)
                call lccree(1, metaRela, metaRelaPy)
                call lccree(1, metaGlob, metaGlobPy)
                call lcinfo(metaPhasPy, idummy, nbMetaPhas, idummy2)
                call lcinfo(metaRelaPy, idummy, nbVariMetaRela, idummy2)
                call lcinfo(metaGlobPy, idummy, nbVariMetaGlob, idummy2)
                ASSERT(nbMetaPhas .le. 10)
                ASSERT(nbVariMetaRela .le. 30)
                ASSERT(nbVariMetaGlob .le. 30)
                call lcvari(metaPhasPy, nbMetaPhas, metaPhasName)
                call lcvari(metaRelaPy, nbVariMetaRela, metaRelaName)
                call lcvari(metaGlobPy, nbVariMetaGlob, metaGlobName)
                iVari = 0
                do iMetaPhas = 1, nbMetaPhas
                    do iVariMetaRela = 1, nbVariMetaRela
                        iVari = iVari + 1
                        infoVari(iVari) = metaPhasName(iMetaPhas)//'##'//metaRelaName(iVariMetaRela)
                    enddo
                enddo
                do iVariMetaRela = 1, nbVariMetaRela
                    iVari = iVari + 1
                    infoVari(iVari) = metaRelaName(iVariMetaRela)
                enddo
                do iVariMetaGlob = 1, nbVariMetaGlob
                    iVari = iVari + 1
                    infoVari(iVari) = metaGlobName(iVariMetaGlob)
                enddo
                call lcdiscard(metaPhasPy)
                call lcdiscard(metaRelaPy)
                call lcdiscard(metaGlobPy)
                nbVariOther = nbVari - iVari
                if (nbVariOther .ne. 0) then
                    call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp,&
                                        post_iter, regu_visc, l_implex,&
                                        compCodePy)
                    call lcvari(compCodePy, nbVariOther, infoVari(iVari+1:nbVari))
                    call lcdiscard(compCodePy)
                endif
            else
                call lcvari(compCodePy, nbVari, infoVari)
                call lcdiscard(compCodePy)
            endif
        endif
    endif
!
end subroutine
