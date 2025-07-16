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
subroutine comp_meca_rkit(keywordfact, iocc, rela_comp, kit_comp, l_etat_init_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/getvtx.h"
#include "asterfort/assert.h"
#include "asterfort/ddi_kit_read.h"
#include "asterfort/thm_kit_read.h"
!
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(out) :: kit_comp(4)
    aster_logical, optional, intent(in) :: l_etat_init_
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Read informations for KIT
!
! --------------------------------------------------------------------------------------------------
!
! In  keywordfact      : factor keyword to read
! In  iocc             : factor keyword index
! In  rela_comp        : comportment relation
! Out kit_comp         : KIT comportment
! In  l_etat_init      : .true. if initial state is defined
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nocc
    character(len=16) :: rela_thmc, rela_hydr, rela_meca, rela_ther
    character(len=16) :: rela_flua, rela_plas, rela_cpla, rela_coup
    character(len=16) :: rela_cg(2)
    character(len=16) :: metaPhas, metaRela, relaCompMeta, metaGlob
    aster_logical :: l_etat_init, lIsot, lCine
!
! --------------------------------------------------------------------------------------------------
!
    kit_comp(1:4) = 'VIDE'
    l_etat_init = .false.
    if (present(l_etat_init_)) then
        l_etat_init = l_etat_init_
    end if
!
    if (rela_comp .eq. 'KIT_META') then
        metaPhas = 'VIDE'
        call getvtx(keywordfact, 'RELATION_KIT', iocc=iocc, &
                    nbval=1, vect=metaPhas, nbret=nocc)
        ASSERT(nocc .eq. 1)

        call getvtx(keywordfact, 'RELATION', iocc=iocc, scal=relaCompMeta)

! ----- Internal state variables (by phase)
        metaRela = 'VIDE'
        lIsot = ASTER_FALSE
        lCine = ASTER_FALSE
        if (relaCompMeta(1:9) .eq. 'META_P_CL') then
            metaRela = 'META_P_CINE_LINE'
            lCine = ASTER_TRUE
        elseif (relaCompMeta(1:9) .eq. 'META_P_IL') then
            metaRela = 'META_P_ISOT_LINE'
            lIsot = ASTER_TRUE
        elseif (relaCompMeta(1:10) .eq. 'META_P_INL') then
            metaRela = 'META_P_ISOT_TRAC'
            lIsot = ASTER_TRUE
        elseif (relaCompMeta(1:9) .eq. 'META_V_CL') then
            metaRela = 'META_V_CINE_LINE'
            lCine = ASTER_TRUE
        elseif (relaCompMeta(1:9) .eq. 'META_V_IL') then
            metaRela = 'META_V_ISOT_LINE'
            lIsot = ASTER_TRUE
        elseif (relaCompMeta(1:10) .eq. 'META_V_INL') then
            metaRela = 'META_V_ISOT_TRAC'
            lIsot = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Internal state variables (global)
        metaGlob = 'VIDE'
        if (lIsot) then
            metaGlob = 'META_G_ISOT'
        elseif (lCine) then
            metaGlob = 'META_G_CINE'
        else
            ASSERT(ASTER_FALSE)
        end if
        if ((relaCompMeta(11:13) .eq. 'PT ') .or. (relaCompMeta(12:14) .eq. 'PT ')) then
            metaGlob(12:16) = '_PT  '
        end if
        if ((relaCompMeta(11:15) .eq. 'PT_RE') .or. (relaCompMeta(12:16) .eq. 'PT_RE')) then
            metaGlob(12:16) = '_PTRE'
        else if ((relaCompMeta(11:13) .eq. 'RE') .or. (relaCompMeta(12:14) .eq. 'RE')) then
            metaGlob(12:16) = '_RE  '
        end if
        kit_comp(1) = metaPhas
        kit_comp(2) = metaRela
        kit_comp(3) = metaGlob

    else if (rela_comp .eq. 'KIT_DDI') then
        call ddi_kit_read(keywordfact, iocc, l_etat_init, &
                          rela_flua, rela_plas, rela_cpla, rela_coup)
        kit_comp(1) = rela_flua
        kit_comp(2) = rela_plas
        kit_comp(3) = rela_coup
        kit_comp(4) = rela_cpla

    else if (rela_comp .eq. 'KIT_CG') then
        call getvtx(keywordfact, 'RELATION_KIT', iocc=iocc, &
                    nbval=2, vect=rela_cg, nbret=nocc)
        ASSERT(nocc .eq. 2)
        kit_comp(1) = rela_cg(1)
        kit_comp(2) = rela_cg(2)

    elseif ((rela_comp(1:5) .eq. 'KIT_H') .or. (rela_comp(1:6) .eq. 'KIT_TH')) then
        call thm_kit_read(keywordfact, iocc, &
                          rela_comp, rela_thmc, rela_hydr, rela_meca, rela_ther)
        kit_comp(1) = rela_meca
        kit_comp(2) = rela_hydr
        kit_comp(3) = rela_ther
        kit_comp(4) = rela_thmc

    else
        ASSERT(ASTER_FALSE)
    end if
end subroutine
