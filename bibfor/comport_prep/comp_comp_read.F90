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

subroutine comp_comp_read(v_info_valk, v_info_vali)
!
!
    implicit none
!
#include "MultiFiber_type.h"
#include "asterc/getfac.h"
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
#include "asterfort/getvid.h"
#include "asterfort/assert.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=16), intent(out) :: v_info_valk(:)
    integer(kind=8), intent(out) :: v_info_vali(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (AFFE_MATERIAU)
!
! Read informations from command file
!
! --------------------------------------------------------------------------------------------------
!
! Out v_info_valk      : comportment informations (character)
! Out v_info_vali      : comportment informations (integer)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: keywordfact
    integer(kind=8) :: iocc, nocc, nbgmax, i, jdecal, idummy, idummy2
    character(len=8) :: sdcomp, k8dummy
    character(len=16) :: comp_code
    character(len=16) :: rela_comp, defo_comp, type_comp, type_cpla, mult_comp
    integer(kind=8) :: nb_vari, nume_comp, nb_vari_exte, unit_comp
    character(len=24), pointer :: v_sdcomp_cprk(:) => null()
    integer(kind=8), pointer :: v_sdcomp_cpri(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    keywordfact = 'AFFE_COMPOR'
    call getfac(keywordfact, nocc)
!
    do iocc = 1, nocc
!
        nb_vari_exte = 0
        unit_comp = 0
!
! ----- Get SD COMPOR from DEFI_COMPOR
!
        call getvid(keywordfact, 'COMPOR', iocc=iocc, scal=sdcomp)
        call jeveuo(sdcomp//'.CPRI', 'L', vi=v_sdcomp_cpri)
        call jeveuo(sdcomp//'.CPRK', 'L', vk24=v_sdcomp_cprk)
        ! Ceinture et bretelle : seulement pour les multifibres
        ASSERT(v_sdcomp_cpri(MULTI_FIBER_TYPE) .eq. 3)
!
! ----- First none-void COMPOR in fiber
!
        call jelira(sdcomp//'.CPRK', 'LONMAX', nbgmax, k8dummy)
        jdecal = 0
        ! -1 car le dernier c'est le matériau de torsion
        nbgmax = (nbgmax-1)/MULTI_FIBER_SIZEK
        do i = 1, nbgmax
            jdecal = MULTI_FIBER_SIZEK*(i-1)
            if (v_sdcomp_cprk(jdecal+MULTI_FIBER_MATER) .ne. 'VIDE') then
                goto 25
            end if
        end do
        call utmess('F', 'COMPOR1_85')
25      continue
!
! ----- Save options in list
!
        rela_comp = v_sdcomp_cprk(jdecal+MULTI_FIBER_RELA) (1:16)
        defo_comp = v_sdcomp_cprk(jdecal+MULTI_FIBER_DEFO) (1:16)
        type_comp = 'COMP_INCR'
        type_cpla = v_sdcomp_cprk(jdecal+MULTI_FIBER_ALGO) (1:16)
        mult_comp = sdcomp//'.CPRK'
        nb_vari = v_sdcomp_cpri(MULTI_FIBER_NBVARMAX)
        call lccree(1, rela_comp, comp_code)
        call lcinfo(comp_code, nume_comp, idummy, idummy2)
        call lcdiscard(comp_code)
!
! ----- Save options in list
!
        v_info_vali(4*(iocc-1)+1) = nb_vari_exte
        v_info_vali(4*(iocc-1)+2) = unit_comp
        v_info_vali(4*(iocc-1)+3) = nb_vari
        v_info_vali(4*(iocc-1)+4) = nume_comp
        v_info_valk(16*(iocc-1)+1) = rela_comp
        v_info_valk(16*(iocc-1)+2) = defo_comp
        v_info_valk(16*(iocc-1)+3) = type_comp
        v_info_valk(16*(iocc-1)+4) = type_cpla
        v_info_valk(16*(iocc-1)+5) = 'VIDE'
        v_info_valk(16*(iocc-1)+6) = 'VIDE'
        v_info_valk(16*(iocc-1)+7) = 'VIDE'
        v_info_valk(16*(iocc-1)+8) = 'VIDE'
        v_info_valk(16*(iocc-1)+14) = mult_comp
        v_info_valk(16*(iocc-1)+15) = 'VIDE'
        v_info_valk(16*(iocc-1)+16) = 'VIDE'
!
    end do
!
end subroutine
