! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
!
subroutine cazocx(sdcont, model, keywf, i_zone)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/exixfe.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: sdcont
    character(len=8), intent(in) :: model
    character(len=16), intent(in) :: keywf
    integer, intent(in) :: i_zone
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Discrete methods - Get parameters of contact zone
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  keywf            : factor keyword to read
! In  i_zone           : index of contact zone
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_defi
    integer :: zcmxf, ztole
    character(len=16) :: s_type_inte, s_algo_lagr, s_gliss, s_czm_rela
    character(len=16) :: s_algo_cont, s_algo_frot, s_cont_init
    integer :: iret, noc, inte_order
    real(kind=8) :: coef_pena_cont, coef_pena_frot, coef_augm_cont, coef_augm_frot
    real(kind=8) :: algo_cont, algo_frot
    real(kind=8) :: coef_coul_frot, seuil_init, tole_proj_ext
    aster_logical :: l_frot
    integer, pointer :: v_xfem_cont(:) => null()
    character(len=24) :: sdcont_caraxf
    real(kind=8), pointer :: v_sdcont_caraxf(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    tole_proj_ext = 0.d0
    coef_pena_cont = 100.d0
    coef_augm_cont = 100.d0
    coef_pena_frot = 100.d0
    coef_augm_frot = 100.d0
    coef_coul_frot = 0.d0
    algo_cont = 0.d0
    algo_frot = 0.d0
    coef_coul_frot = 0.d0
    seuil_init = 0.d0
    s_type_inte = 'FPG4'
    s_algo_lagr = 'NON'
    s_algo_cont = 'STANDARD'
    s_algo_frot = 'STANDARD'
    s_gliss = ' '
    s_czm_rela = ' '
    s_cont_init = ' '
!
! - Datastructure for contact definition
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_caraxf = sdcont_defi(1:16)//'.CARAXF'
    call jeveuo(sdcont_caraxf, 'E', vr=v_sdcont_caraxf)
    zcmxf = cfmmvd('ZCMXF')
    ztole = cfmmvd('ZTOLE')
!
! - Parameters
!
    l_frot = cfdisl(sdcont_defi, 'FROTTEMENT')
!
! - Check model (XFEM)
!
    call exixfe(model, iret)
    if (iret .eq. 0) then
        call utmess('F', 'XFEM2_8')
    else
        call jeveuo(model(1:8)//'.XFEM_CONT', 'L', vi=v_xfem_cont)
        if (v_xfem_cont(1) .eq. 0) then
            call utmess('F', 'XFEM2_9')
        end if
    end if
!
! - Integration scheme
!
    call getvtx(keywf, 'INTEGRATION', iocc=i_zone, scal=s_type_inte)
    if (s_type_inte .eq. 'NOEUD') then
        v_sdcont_caraxf(zcmxf*(i_zone-1)+1) = 1.d0
    else if (s_type_inte .eq. 'GAUSS') then
        call getvis(keywf, 'ORDRE_INT', iocc=i_zone, scal=inte_order)
        v_sdcont_caraxf(zcmxf*(i_zone-1)+1) = 10.d0*inte_order+2.d0
    else if (s_type_inte(1:7) .eq. 'SIMPSON') then
        call getvis(keywf, 'ORDRE_INT', iocc=i_zone, scal=inte_order)
        v_sdcont_caraxf(zcmxf*(i_zone-1)+1) = 10.d0*inte_order+3.d0
    else if (s_type_inte(1:6) .eq. 'NCOTES') then
        call getvis(keywf, 'ORDRE_INT', iocc=i_zone, scal=inte_order)
        v_sdcont_caraxf(zcmxf*(i_zone-1)+1) = 10.d0*inte_order+4.d0
    else
        ASSERT(.false.)
    end if
!
! - Contact method
!
    call getvtx(keywf, 'ALGO_CONT', iocc=i_zone, scal=s_algo_cont)
    if (s_algo_cont .eq. 'STANDARD') then
        call getvr8(keywf, 'COEF_CONT', iocc=i_zone, scal=coef_augm_cont)
        algo_cont = 1.d0
        coef_pena_cont = 0.d0
    else if (s_algo_cont .eq. 'PENALISATION') then
        call getvr8(keywf, 'COEF_PENA_CONT', iocc=i_zone, scal=coef_pena_cont)
        algo_cont = 2.d0
        coef_augm_cont = 0.d0
    else if (s_algo_cont .eq. 'CZM') then
        algo_cont = 3.d0
        l_frot = .false.
    else
        ASSERT(.false.)
    end if
    v_sdcont_caraxf(zcmxf*(i_zone-1)+2) = coef_augm_cont
    v_sdcont_caraxf(zcmxf*(i_zone-1)+12) = coef_pena_cont
    v_sdcont_caraxf(zcmxf*(i_zone-1)+11) = algo_cont
!
! - Friction method
!
    if (l_frot) then
        call getvtx(keywf, 'ALGO_FROT', iocc=i_zone, scal=s_algo_frot)
        if (s_algo_frot .eq. 'STANDARD') then
            call getvr8(keywf, 'COEF_FROT', iocc=i_zone, scal=coef_augm_frot)
            algo_frot = 1.d0
            coef_pena_frot = 0.d0
        else if (s_algo_frot(1:14) .eq. 'PENALISATION') then
            call getvr8(keywf, 'COEF_PENA_FROT', iocc=i_zone, scal=coef_pena_frot)
            algo_frot = 2.d0
            coef_augm_frot = 0.d0
        else
            ASSERT(.false.)
        end if
    else
        coef_augm_frot = 0.d0
        coef_pena_frot = 0.d0
        algo_frot = 0.d0
    end if
    v_sdcont_caraxf(zcmxf*(i_zone-1)+3) = coef_augm_frot
    v_sdcont_caraxf(zcmxf*(i_zone-1)+14) = coef_pena_frot
    v_sdcont_caraxf(zcmxf*(i_zone-1)+13) = algo_frot
!
! - Parameters of friction
!
    if (l_frot) then
        v_sdcont_caraxf(zcmxf*(i_zone-1)+5) = 3.d0
        call getvr8(keywf, 'COULOMB', iocc=i_zone, scal=coef_coul_frot)
        v_sdcont_caraxf(zcmxf*(i_zone-1)+4) = coef_coul_frot
        call getvr8(keywf, 'SEUIL_INIT', iocc=i_zone, scal=seuil_init)
        v_sdcont_caraxf(zcmxf*(i_zone-1)+6) = seuil_init
    else
        v_sdcont_caraxf(zcmxf*(i_zone-1)+5) = 1.d0
    end if
!
! - Parameters of contact
!
    if (s_algo_cont .ne. 'CZM') then
!
! ----- Initial contact
!
        call getvtx(keywf, 'CONTACT_INIT', iocc=i_zone, scal=s_cont_init, nbret=noc)
        if (s_cont_init .eq. 'OUI') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+7) = 1.d0
        else if (s_cont_init .eq. 'NON') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+7) = 0.d0
        else
            ASSERT(.false.)
        end if
!
! ----- Bilateral contact
!
        call getvtx(keywf, 'GLISSIERE', iocc=i_zone, scal=s_gliss)
        if (s_gliss .eq. 'OUI') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+10) = 1.d0
        else if (s_gliss .eq. 'NON') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+10) = 0.d0
        else
            ASSERT(.false.)
        end if
    end if
!
! - Parameters of Lagrange multipliers algorithm
!
    call getvtx(keywf, 'ALGO_LAGR', iocc=1, scal=s_algo_lagr)
    if (s_algo_lagr .eq. 'NON') then
        v_sdcont_caraxf(zcmxf*(i_zone-1)+9) = 0.d0
    else if (s_algo_lagr .eq. 'VERSION1') then
        v_sdcont_caraxf(zcmxf*(i_zone-1)+9) = 1.d0
    else if (s_algo_lagr .eq. 'VERSION2') then
        v_sdcont_caraxf(zcmxf*(i_zone-1)+9) = 2.d0
    else if (s_algo_lagr .eq. 'VERSION3') then
!       Algorithme uniquement disponible pour des éléments quadratiques
        if (v_xfem_cont(1) .eq. 1) call utmess('F', 'XFEM_90')

        v_sdcont_caraxf(zcmxf*(i_zone-1)+9) = 3.d0
    else if (s_algo_lagr .eq. 'AUTO') then
!       cas 'AUTO' :
!          - si quadratique et formulation petits glissements,
!            on choisit l'algo P2/P1*
!          - sinon, on choisit l'algo "version 2" (P1/P1 ou P2/P1-)
        if (v_xfem_cont(1) .eq. 3) then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+9) = 3.d0
        else
            v_sdcont_caraxf(zcmxf*(i_zone-1)+9) = 2.d0
        end if
    else
        ASSERT(.false.)
    end if
!
! - TOLE_PROJ_EXT
! --- TOLE_PROJ_EXT <0: disallow projection outside element
! --- TOLE_PROJ_EXT >0: allow projection outside element with TOLE_PROJ_EXT tolerance
!
    call getvr8(keywf, 'TOLE_PROJ_EXT', iocc=i_zone, scal=tole_proj_ext)
    if (tole_proj_ext .lt. 0.d0) then
        v_sdcont_caraxf(zcmxf*(i_zone-1)+15) = -1.d0
    else
        v_sdcont_caraxf(zcmxf*(i_zone-1)+15) = tole_proj_ext
    end if
!
! - CZM parameters
!
    if (s_algo_cont .eq. 'CZM') then
        call getvtx(keywf, 'RELATION', iocc=i_zone, scal=s_czm_rela)
        if (s_czm_rela .eq. 'CZM_EXP_REG') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+16) = 1.d0
        else if (s_czm_rela .eq. 'CZM_LIN_REG') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+16) = 2.d0
        else if (s_czm_rela .eq. 'CZM_TAC_MIX') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+16) = 3.d0
        else if (s_czm_rela .eq. 'CZM_OUV_MIX') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+16) = 4.d0
        else if (s_czm_rela .eq. 'CZM_LIN_MIX') then
            v_sdcont_caraxf(zcmxf*(i_zone-1)+16) = 5.d0
        else
            v_sdcont_caraxf(zcmxf*(i_zone-1)+16) = 0.d0
        end if
!
! ----- Check
!
        if (s_czm_rela .eq. 'CZM_LIN_MIX' .and. v_xfem_cont(1) .ne. 2) then
            call utmess('F', 'XFEM_93', sk=model)
        else if (s_czm_rela .ne. 'CZM_LIN_MIX' .and. v_xfem_cont(1) .ne. 1 .and. &
                 v_xfem_cont(1) .ne. 3) then
            call utmess('F', 'XFEM_93', sk=model)
        end if
    end if
!
end subroutine
