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
subroutine comp_meca_name(nbVari, nbVariMeca, l_excl, vari_excl, l_kit_meta, &
                          rela_comp, defo_comp, kit_comp, type_cpla, post_iter, &
                          regu_visc, post_incr, &
                          extern_addr, extern_type, infoVari)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcvari.h"
#include "asterc/lcdiscard.h"
#include "asterfort/assert.h"
#include "asterfort/comp_mfront_vname.h"
#include "asterfort/comp_meca_code.h"
!
    integer(kind=8), intent(in) :: nbVari, nbVariMeca
    aster_logical, intent(in) :: l_excl
    character(len=16), intent(in) :: vari_excl
    aster_logical, intent(in) :: l_kit_meta
    character(len=16), intent(in) :: extern_addr, rela_comp, defo_comp, kit_comp(4)
    character(len=16), intent(in) :: type_cpla, post_iter, regu_visc, post_incr
    integer(kind=8), intent(in) :: extern_type
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
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  kit_comp         : KIT comportment
! In  type_cpla        : plane stress method
! In  post_iter        : type of post-treatment at each Newton iteration
! In  regu_visc        : keyword for viscuous regularization
! In  post_incr        : type of post-treatment at end of time step
! Ptr infoVari         : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6) :: metaPhasName(10)
    character(len=8) :: metaRelaName(30)
    character(len=16) :: metaGlobName(30)
    integer(kind=8) :: idummy, idummy2, nbVariOther, iVariMeca, iVari
    character(len=16) :: compCodePy
    character(len=16) :: metaPhas, metaRela, metaGlob
    character(len=16) :: metaPhasPy, metaRelaPy, metaGlobPy
    integer(kind=8) :: nbMetaPhas, nbVariMetaRela, nbVariMetaGlob
    integer(kind=8) :: iMetaPhas, iVariMetaRela, iVariMetaGlob
!
! --------------------------------------------------------------------------------------------------
!
    if (l_excl) then
        infoVari(1:nbVari) = vari_excl
    else
! ----- Name of internal state variables
        ! UMAT
        if (extern_type .eq. 4) then
            call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp, &
                                post_iter, regu_visc, post_incr, &
                                compCodePy)
            nbVariOther = nbVari-nbVariMeca
            do iVariMeca = 1, nbVariMeca
                infoVari(iVariMeca) = 'NoName'
            end do
            if (nbVariOther .ne. 0) then
                call lcvari(compCodePy, nbVariOther, infoVari(nbVariMeca+1:nbVari))
            end if
            call lcdiscard(compCodePy)

            ! MFront official or proto
        else if (extern_type .eq. 1 .or. extern_type .eq. 2) then
            ASSERT(extern_addr .ne. ' ')
            call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp, &
                                post_iter, regu_visc, post_incr, &
                                compCodePy)
            nbVariOther = nbVari-nbVariMeca
            call comp_mfront_vname(extern_addr, nbVariMeca, infoVari)
            if (nbVariOther .ne. 0) then
                call lcvari(compCodePy, nbVariOther, infoVari(nbVariMeca+1:nbVari))
            end if
            call lcdiscard(compCodePy)

            ! internal integration
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
                        iVari = iVari+1
                        infoVari(iVari) = metaPhasName(iMetaPhas)//'##'//metaRelaName(iVariMetaRela)
                    end do
                end do
                do iVariMetaGlob = 1, nbVariMetaGlob
                    iVari = iVari+1
                    infoVari(iVari) = metaGlobName(iVariMetaGlob)
                end do
                ASSERT(iVari .eq. nbVariMetaGlob+nbVariMetaRela*nbMetaPhas)
                call lcdiscard(metaPhasPy)
                call lcdiscard(metaRelaPy)
                call lcdiscard(metaGlobPy)
                nbVariOther = nbVari-iVari
                if (nbVariOther .ne. 0) then
                    call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp, &
                                        post_iter, regu_visc, post_incr, &
                                        compCodePy)
                    call lcvari(compCodePy, nbVariOther, infoVari(iVari+1:nbVari))
                    call lcdiscard(compCodePy)
                end if
            else
                call comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp, &
                                    post_iter, regu_visc, post_incr, &
                                    compCodePy)
                call lcvari(compCodePy, nbVari, infoVari)
                call lcdiscard(compCodePy)
            end if
        end if
    end if
!
end subroutine
