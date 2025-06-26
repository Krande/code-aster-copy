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
!
subroutine cazocp(sdcont)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "Contact_type.h"
!
    character(len=8), intent(in) :: sdcont
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get parameters (not depending on contact zones)
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont         : name of contact datastructure
! In  model          : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_defi
    integer(kind=8) :: geom_nbiter, nb_resol, gcp_maxi, gcp_pre_maxi, cont_stat_elas
    integer(kind=8) :: cont_mult, frot_maxi, geom_maxi
    character(len=16) :: gcp_rech_line, gcp_precond, reac_geom, cont_type, stop_singular, elim_edge
    character(len=16) :: algo_reso_cont, algo_reso_frot, algo_reso_geom
    integer(kind=8) :: noc
    real(kind=8) :: resi_abso, gcp_coef_resi
    real(kind=8) :: geom_resi, frot_resi, cont_resi
    aster_logical :: l_cont_gcp
    aster_logical :: l_cont_disc, l_cont_cont, l_frot, l_cont_lac
    character(len=16) :: lissage
    character(len=24) :: sdcont_paracr
    real(kind=8), pointer :: v_sdcont_paracr(:) => null()
    character(len=24) :: sdcont_paraci
    integer(kind=8), pointer :: v_sdcont_paraci(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    elim_edge = 'DUAL'
    reac_geom = 'AUTOMATIQUE'
    algo_reso_cont = ' '
    algo_reso_frot = ' '
    algo_reso_geom = ' '
    geom_nbiter = 2
    nb_resol = 10
    geom_resi = 1.d-2
    frot_resi = 1.d-2
    cont_resi = 1.d-2
!
! - Access to datastructure
!
    sdcont_paracr = sdcont(1:8)//'.PARACR'
    sdcont_paraci = sdcont(1:8)//'.PARACI'
    call jeveuo(sdcont_paracr, 'E', vr=v_sdcont_paracr)
    call jeveuo(sdcont_paraci, 'E', vi=v_sdcont_paraci)
!
! - Active functionnalites
!
    l_cont_disc = cfdisl(sdcont_defi, 'FORMUL_DISCRETE')
    l_cont_cont = cfdisl(sdcont_defi, 'FORMUL_CONTINUE')
    l_cont_lac = cfdisl(sdcont_defi, 'FORMUL_LAC')
    l_cont_gcp = cfdisl(sdcont_defi, 'CONT_GCP')
    l_frot = cfdisl(sdcont_defi, 'FROTTEMENT')
!
! - Geometric algorithm
!
    if (l_cont_cont) then
        call getvtx(' ', 'ALGO_RESO_GEOM', scal=algo_reso_geom)
    else if (l_cont_disc) then
        algo_reso_geom = 'POINT_FIXE'
    else if (l_cont_lac) then
        call getvtx(' ', 'ALGO_RESO_GEOM', scal=algo_reso_geom)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (algo_reso_geom .eq. 'POINT_FIXE') then
        v_sdcont_paraci(9) = ALGO_FIXE
    else if (algo_reso_geom .eq. 'NEWTON') then
        v_sdcont_paraci(9) = ALGO_NEWT
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Geometric parameters
!
    if (algo_reso_geom .eq. 'POINT_FIXE') then
        call getvtx(' ', 'REAC_GEOM', scal=reac_geom)
        if (reac_geom .eq. 'SANS') then
            v_sdcont_paraci(1) = 0
            v_sdcont_paracr(1) = geom_resi
        else if (reac_geom .eq. 'AUTOMATIQUE') then
            v_sdcont_paraci(1) = -1
            call getvis(' ', 'ITER_GEOM_MAXI', scal=geom_maxi)
            v_sdcont_paraci(6) = geom_maxi
            call getvr8(' ', 'RESI_GEOM', scal=geom_resi)
            v_sdcont_paracr(1) = geom_resi
        else if (reac_geom .eq. 'CONTROLE') then
            call getvis(' ', 'NB_ITER_GEOM', scal=geom_nbiter)
            v_sdcont_paraci(1) = geom_nbiter
            v_sdcont_paracr(1) = geom_resi
        else
            ASSERT(ASTER_FALSE)
        end if
        if (l_cont_lac) then
            call utmess('F', 'CONTACT4_1')
        end if
    else if (algo_reso_geom .eq. 'NEWTON') then
        call getvr8(' ', 'RESI_GEOM', scal=geom_resi)
        v_sdcont_paraci(1) = 0
        v_sdcont_paracr(1) = geom_resi
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Friction algorithm
!
    if (l_frot) then
        if (l_cont_cont) then
            call getvtx(' ', 'ALGO_RESO_FROT', scal=algo_reso_frot)
        else if (l_cont_disc) then
            algo_reso_frot = 'POINT_FIXE'
        else if (l_cont_lac) then
            call utmess('F', 'CONTACT4_4')
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    if (l_frot) then
        if (algo_reso_frot .eq. 'POINT_FIXE') then
            v_sdcont_paraci(28) = ALGO_FIXE
        else if (algo_reso_frot .eq. 'NEWTON') then
            v_sdcont_paraci(28) = ALGO_NEWT
        else
            ASSERT(ASTER_FALSE)
        end if
        ASSERT(.not. l_cont_lac)
    end if
!
! - Friction parameters
!
    if (l_frot) then
        if (l_cont_cont) then
            if (algo_reso_frot .eq. 'POINT_FIXE') then
                call getvis(' ', 'ITER_FROT_MAXI', scal=frot_maxi)
                v_sdcont_paraci(7) = frot_maxi
                call getvr8(' ', 'RESI_FROT', scal=frot_resi)
                v_sdcont_paracr(2) = frot_resi
            else
                call getvr8(' ', 'RESI_FROT', scal=frot_resi)
                v_sdcont_paracr(2) = frot_resi
            end if
        end if
        ASSERT(.not. l_cont_lac)
    end if
!
! - Contact algorithm
!
    if (l_cont_cont) then
        call getvtx(' ', 'ALGO_RESO_CONT', scal=algo_reso_cont)
    else if (l_cont_disc) then
        algo_reso_cont = 'POINT_FIXE'
    else if (l_cont_lac) then
        call getvtx(' ', 'ALGO_RESO_CONT', scal=algo_reso_cont)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (algo_reso_cont .eq. 'POINT_FIXE') then
        v_sdcont_paraci(27) = ALGO_FIXE
    else if (algo_reso_cont .eq. 'NEWTON') then
        v_sdcont_paraci(27) = ALGO_NEWT
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Contact status stabilisation with MATRICE="ELASTQUE"
!
    if (l_cont_cont .or. l_cont_lac) then
        call getvis(' ', 'CONT_STAT_ELAS', scal=cont_stat_elas)
        v_sdcont_paraci(31) = cont_stat_elas
    end if
!
! - Check
!
    if (algo_reso_geom .eq. 'NEWTON' .and. algo_reso_cont .ne. 'NEWTON') then
        call utmess('F', 'CONTACT_20')
    end if
!
! - Contact parameters
!
    if (algo_reso_cont .eq. 'POINT_FIXE') then
        if (l_cont_cont) then
            call getvtx(' ', 'ITER_CONT_TYPE', scal=cont_type)
            if (cont_type .eq. 'MULT') then
                cont_mult = 4
                call getvis(' ', 'ITER_CONT_MULT', scal=cont_mult)
                v_sdcont_paraci(5) = cont_mult
                v_sdcont_paraci(10) = -1
            else if (cont_type .eq. 'MAXI') then
                cont_mult = 30
                call getvis(' ', 'ITER_CONT_MAXI', scal=cont_mult)
                v_sdcont_paraci(10) = cont_mult
                v_sdcont_paraci(5) = -1
            else
                ASSERT(ASTER_FALSE)
            end if
            v_sdcont_paracr(7) = -1
        else if (l_cont_disc) then
            call getvis(' ', 'ITER_CONT_MULT', scal=cont_mult)
            v_sdcont_paraci(5) = cont_mult
            v_sdcont_paraci(10) = -1
        else if (l_cont_lac) then
            call utmess('F', 'CONTACT4_1')
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (algo_reso_cont .eq. 'NEWTON') then
        if (l_cont_cont) then
            call getvr8(' ', 'RESI_CONT', scal=cont_resi)
            v_sdcont_paracr(7) = cont_resi
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Discrete formulation
!
    if (l_cont_disc) then
        call getvtx(' ', 'STOP_SINGULIER', scal=stop_singular)
        if (stop_singular .eq. 'OUI') then
            v_sdcont_paraci(2) = 0
        else if (stop_singular .eq. 'NON') then
            v_sdcont_paraci(2) = 1
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call getvis(' ', 'NB_RESOL', scal=nb_resol)
        v_sdcont_paraci(3) = nb_resol
!
        if (l_cont_gcp) then
            call getvr8(' ', 'RESI_ABSO', scal=resi_abso, nbret=noc)
            if (noc .eq. 0) then
                call utmess('F', 'CONTACT_4')
            end if
            v_sdcont_paracr(4) = resi_abso
!
            call getvis(' ', 'ITER_GCP_MAXI', scal=gcp_maxi)
            v_sdcont_paraci(12) = gcp_maxi
!
            call getvtx(' ', 'PRE_COND', scal=gcp_precond)
            if (gcp_precond .eq. 'SANS') then
                v_sdcont_paraci(13) = 0
            else if (gcp_precond .eq. 'DIRICHLET') then
                v_sdcont_paraci(13) = 1
                call getvr8(' ', 'COEF_RESI', scal=gcp_coef_resi)
                v_sdcont_paracr(5) = gcp_coef_resi
                call getvis(' ', 'ITER_PRE_MAXI', scal=gcp_pre_maxi)
                v_sdcont_paraci(14) = gcp_pre_maxi
            else
                ASSERT(ASTER_FALSE)
            end if
!
            call getvtx(' ', 'RECH_LINEAIRE', scal=gcp_rech_line)
            if (gcp_rech_line .eq. 'ADMISSIBLE') then
                v_sdcont_paraci(15) = 0
            else if (gcp_rech_line .eq. 'NON_ADMISSIBLE') then
                v_sdcont_paraci(15) = 1
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    end if
!
! - Smoothing
!
    if (l_cont_disc .or. l_cont_cont .or. l_cont_lac) then
        call getvtx(' ', 'LISSAGE', scal=lissage)
        if (lissage(1:3) .eq. 'NON') then
            v_sdcont_paraci(19) = 0
        else if (lissage(1:3) .eq. 'OUI') then
            v_sdcont_paraci(19) = 1
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Verification method
    if (l_cont_disc .or. l_cont_cont) then
        call getvtx(' ', 'STOP_INTERP', scal=stop_singular)
        if (stop_singular .eq. 'OUI') then
            v_sdcont_paraci(25) = 1
        else if (stop_singular .eq. 'NON') then
            v_sdcont_paraci(25) = 0
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
end subroutine
