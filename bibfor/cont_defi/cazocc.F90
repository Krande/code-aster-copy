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
subroutine cazocc(sdcont, factorKeyword, i_zone, nb_cont_zone)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfl.h"
#include "asterfort/utmess.h"
#include "asterc/r8prem.h"
!
    character(len=8), intent(in) :: sdcont
    integer(kind=8), intent(in) :: i_zone
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), optional, intent(in) :: nb_cont_zone
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Continue method - Get parameters of contact zone
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  factorKeyword            : factor keyword to read
! In  i_zone           : index of contact zone
! In  nb_cont_zone     : number of zones of contact
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_defi
    integer(kind=8) :: nb_frot_excl_1, nb_frot_excl_2
    integer(kind=8) :: nb_cont_excl_1, nb_cont_excl_2, nb_cont_excl_3, nb_cont_excl_4
    integer(kind=8) :: nb_dire_excl, noc, iadapt
    character(len=16) :: s_cont_excl, s_frot_excl
    character(len=16) :: s_grglis
    character(len=16) :: s_gliss, s_type_inte, s_cont_init, s_algo_cont, s_algo_frot
    character(len=16) :: adaptation
    real(kind=8) :: dire_excl_frot_i, dire_excl_frot(3)
    real(kind=8) :: coef_cont, coef_frot, seuil_init, coef_coul_frot
    real(kind=8) :: coef_augm_frot, coef_augm_cont
    real(kind=8) :: coef_pena_frot, coef_pena_cont, pene_maxi, glis_maxi
    real(kind=8) :: algo_cont, algo_frot
    real(kind=8) :: type_inte, contInit, seuil_auto, contInitDist
    integer(kind=8) :: inte_order
    integer(kind=8) :: nbret
    aster_logical :: l_inte_node, l_frot, l_node_excl, l_frot_excl, l_dire_excl_frot
    aster_logical :: l_gliss
    integer(kind=8) :: zcmcf, zexcl
    character(len=24) :: sdcont_caracf
    real(kind=8), pointer :: v_sdcont_caracf(:) => null()
    character(len=24) :: sdcont_exclfr
    real(kind=8), pointer :: v_sdcont_exclfr(:) => null()
    character(len=24) :: sdcont_paraci
    integer(kind=8), pointer :: v_sdcont_paraci(:) => null()
    aster_logical ::  l_newt_fr, l_cont_cont
    aster_logical ::  l_granglis
    character(len=24) :: sdcont_paracr
    real(kind=8), pointer :: v_sdcont_paracr(:) => null()
    real(kind=8) :: pene_critere
    integer(kind=8) :: i_pene_zone
!
! --------------------------------------------------------------------------------------------------
!
    inte_order = 0
    type_inte = 0.d0
    algo_cont = 0.d0
    coef_augm_cont = 100.d0
    coef_pena_cont = 100.d0
    coef_cont = 100.d0
    algo_frot = 0.d0
    coef_coul_frot = 0.d0
    coef_augm_frot = 100.d0
    coef_pena_frot = 100.d0
    coef_frot = 0.d0
    coef_coul_frot = 0.d0
    seuil_init = 0.d0
    seuil_auto = 0.d0
    contInit = 0.d0
    dire_excl_frot_i = 0.d0
    dire_excl_frot(:) = 0.d0
    l_inte_node = ASTER_FALSE
    l_node_excl = ASTER_FALSE
    l_frot_excl = ASTER_FALSE
    l_gliss = ASTER_FALSE
    l_dire_excl_frot = ASTER_FALSE
    l_newt_fr = ASTER_FALSE
    l_cont_cont = ASTER_FALSE
    l_granglis = ASTER_FALSE
    s_algo_cont = ' '
    s_algo_frot = ' '
    iadapt = 0
    pene_maxi = 1.d3
    glis_maxi = 1.d3
    nbret = -3
    contInitDist = -1.d0

! - Datastructure for contact
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_caracf = sdcont_defi(1:16)//'.CARACF'
    sdcont_exclfr = sdcont_defi(1:16)//'.EXCLFR'
    sdcont_paraci = sdcont(1:8)//'.PARACI'
    sdcont_paracr = sdcont(1:8)//'.PARACR'
    call jeveuo(sdcont_caracf, 'E', vr=v_sdcont_caracf)
    call jeveuo(sdcont_exclfr, 'E', vr=v_sdcont_exclfr)
    call jeveuo(sdcont_paraci, 'E', vi=v_sdcont_paraci)
    call jeveuo(sdcont_paracr, 'E', vr=v_sdcont_paracr)
    zcmcf = cfmmvd('ZCMCF')
    zexcl = cfmmvd('ZEXCL')

! - Parameters
    l_frot = cfdisl(sdcont_defi, 'FROTTEMENT')
    l_cont_cont = cfdisl(sdcont_defi, 'FORMUL_CONTINUE')
    l_newt_fr = cfdisl(sdcont_defi, 'FROT_NEWTON')
!
! - Integration scheme
!
    call getvtx(factorKeyword, 'INTEGRATION', iocc=i_zone, scal=s_type_inte)
    if (s_type_inte .eq. 'AUTO') then
        l_inte_node = ASTER_TRUE
        type_inte = 1.d0
    else if (s_type_inte .eq. 'GAUSS') then
        call getvis(factorKeyword, 'ORDRE_INT', iocc=i_zone, scal=inte_order)
        type_inte = 10.d0*inte_order+2.d0
    else if (s_type_inte .eq. 'SIMPSON') then
        call getvis(factorKeyword, 'ORDRE_INT', iocc=i_zone, scal=inte_order)
        type_inte = 10.d0*inte_order+3.d0
    else if (s_type_inte .eq. 'NCOTES') then
        call getvis(factorKeyword, 'ORDRE_INT', iocc=i_zone, scal=inte_order)
        type_inte = 10.d0*inte_order+4.d0
    else
        ASSERT(ASTER_FALSE)
    end if
    v_sdcont_caracf(zcmcf*(i_zone-1)+1) = type_inte
!
! - Get algorithm for contact
!
    call getvtx(factorKeyword, 'ALGO_CONT', iocc=i_zone, scal=s_algo_cont)
    if (s_algo_cont .eq. 'STANDARD') then
        algo_cont = 1.d0
    else if (s_algo_cont .eq. 'PENALISATION') then
        algo_cont = 3.d0
    else if (s_algo_cont .eq. 'LAC') then
        algo_cont = 5.d0
    else
        ASSERT(ASTER_FALSE)
    end if
    v_sdcont_caracf(zcmcf*(i_zone-1)+3) = algo_cont
!
! - Get algorithm for friction
!
    if (l_frot) then
        call getvtx(factorKeyword, 'ALGO_FROT', iocc=i_zone, scal=s_algo_frot)
        if (s_algo_frot .eq. 'STANDARD') then
            algo_frot = 1.d0
        else if (s_algo_frot .eq. 'PENALISATION') then
            algo_frot = 3.d0
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        algo_frot = 0.d0
    end if
!
! - Automatic adaptation method
!
    call getvtx(factorKeyword, 'ADAPTATION', iocc=i_zone, scal=adaptation, nbret=nbret)
    if (nbret .le. 0) then
        adaptation = 'NON'
    end if
    if (adaptation .eq. 'NON') then
        iadapt = 0
    else if (adaptation .eq. 'ADAPT_COEF') then
        iadapt = -1
! IL FAUT DISTINGUER 3 CAS DE FIGURES :
! FROTTEMENT NEWTON ACTIF : ON ADAPTE COEF_FROT COMME DECRIT DANS LA DOC R + THESE DK
!                    PARACI = 1
!                    SI PENALISATION CONTACT COEF_CONT ADAPTE SELON BUSSETTA
!                    PARACI = 2
! CONTACT PENALISE ACTIF FROTTEMENT TRESCA/NEWTON OU SANS : COEF_CONT ADAPTE SELON BUSSETTA
!                    PARACI = 3
! CONTACT STANDARD ACTIF FROTTEMENT TRESCA ACTIF  : ON NE FAIT RIEN
!                    PARACI = 0
        if (l_cont_cont) then
            if (l_newt_fr) then
                iadapt = 1
                if (s_algo_cont .eq. 'PENALISATION') then
                    iadapt = 2
                end if
                if (s_algo_frot .eq. 'PENALISATION' .and. s_algo_cont .eq. 'STANDARD') then
                    iadapt = 7
                end if
            elseif (s_algo_cont .eq. 'PENALISATION') then
                iadapt = 3
            else
                iadapt = 0
            end if
        end if
    else if (adaptation .eq. 'CYCLAGE') then
        iadapt = 4
    else if (adaptation .eq. 'TOUT') then
        ! IL FAUT DISTINGUER 3 CAS DE FIGURES :
        ! FROTTEMENT NEWTON ACTIF : ON ADAPTE COEF_FROT COMME DECRIT DANS LA DOC R + THESE DK
        !                    PARACI = 1+4
        !                    SI PENALISATION CONTACT COEF_CONT ADAPTE SELON BUSSETTA
        !                    PARACI = 2+4
        ! CONTACT PENALISE ACTIF  FROTTEMENT TRESCA OU NON   : COEF_CONT ADAPTE SELON BUSSETTA
        !                    PARACI = 3+4
        ! CONTACT STANDARD ACTIF FROTTEMENT TRESCA OU NON  : ON NE FAIT RIEN
        !                    PARACI = 0+4
        iadapt = 4
        if (l_cont_cont) then
            if (l_newt_fr) then
                iadapt = iadapt+1
                if (s_algo_cont .eq. 'PENALISATION') then
                    iadapt = iadapt+2
                end if
                if (s_algo_frot .eq. 'PENALISATION' .and. s_algo_cont .eq. 'STANDARD') then
                    iadapt = iadapt+7
                end if
            else if (s_algo_cont .eq. 'PENALISATION') then
                iadapt = iadapt+3
            else
                iadapt = iadapt+0
            end if
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
    v_sdcont_paraci(20) = iadapt
!
! - Contact method
! - Traitement de PENE_MAXI :
!      - Dans le cas  ADAPTATION=NON ou CYCLAGE+PENALISATION,
!        le critere PENE_MAXI n'a pas de sens. Au moment de la résolution,
!        dans STAT_NON_LINE, on ne cherche pas à le vérifier. La parade c'est
!        de prendre le critère volontairement tres grand pour ne pas avoir à le vérifier.
!      - Dans le cas STANDARD,LAC PENE_MAXI n'a pas de sens. Il ne sert qu'à initialiser
!        la variable .PARACR(6) sinon bug dans certaines configurations multi-zone.
!
    if (s_algo_cont .eq. 'STANDARD') then
        call getvr8(factorKeyword, 'COEF_CONT', iocc=i_zone, scal=coef_augm_cont)
        coef_cont = coef_augm_cont
        pene_maxi = 1.d3
    else if (s_algo_cont .eq. 'PENALISATION') then
        ! L'utilisateur peut ne pas renseigner pene_maxi
        call getvr8(factorKeyword, 'PENE_MAXI', iocc=i_zone, scal=pene_maxi, nbret=nbret)
        if ((adaptation .eq. 'NON' .or. adaptation .eq. 'CYCLAGE') .and. (nbret .le. 0)) then
            call getvr8(factorKeyword, 'COEF_PENA_CONT', iocc=i_zone, scal=coef_pena_cont)
            pene_maxi = 1.d3
        elseif (adaptation .eq. 'ADAPT_COEF' .or. adaptation .eq. 'TOUT' .or. (nbret .ge. 1)) then
            if (nbret .le. 0) then
                pene_maxi = -1
            end if
            coef_pena_cont = 1.d2
        end if
        coef_cont = coef_pena_cont
    else if (s_algo_cont .eq. 'LAC') then
        coef_cont = coef_augm_cont
        pene_maxi = 1.d3
        if (iadapt .ne. 0) then
            call utmess('F', 'CONTACT_95')
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
    v_sdcont_caracf(zcmcf*(i_zone-1)+2) = coef_cont
    v_sdcont_caracf(zcmcf*(i_zone-1)+14) = pene_maxi
!
! - Set PENE_MAXI as mininum on all zones
!
    if (s_algo_cont .ne. 'LAC') then
        if (i_zone .eq. nb_cont_zone) then
            pene_critere = pene_maxi
            do i_pene_zone = 1, nb_cont_zone
                if (v_sdcont_caracf(zcmcf*(i_pene_zone-1)+14) .lt. pene_critere) then
                    v_sdcont_paracr(6) = v_sdcont_caracf(zcmcf*(i_pene_zone-1)+14)
                end if
            end do
        end if
    end if
!
! - Get friction parameters
!
    coef_frot = 0.d0
    if (l_frot) then
        if (s_algo_frot .eq. 'STANDARD') then
            call getvr8(factorKeyword, 'COEF_FROT', iocc=i_zone, scal=coef_augm_frot)
            coef_frot = coef_augm_frot
        else if (s_algo_frot .eq. 'PENALISATION') then
            call getvr8(factorKeyword, 'COEF_PENA_FROT', iocc=i_zone, scal=coef_pena_frot)
            call getvr8(factorKeyword, 'GLIS_MAXI', iocc=i_zone, scal=glis_maxi, nbret=nbret)
            coef_frot = coef_pena_frot
        else
            ASSERT(ASTER_FALSE)
        end if
        call getvr8(factorKeyword, 'COULOMB', iocc=i_zone, scal=coef_coul_frot)
        call getvr8(factorKeyword, 'SEUIL_INIT', iocc=i_zone, scal=seuil_init, nbret=noc)
        if (noc .eq. 0) then
            seuil_auto = 1.d0
        end if
        if (coef_coul_frot .le. r8prem()) then
            coef_frot = 0.d0
            algo_frot = 0.d0
        end if
    end if
    v_sdcont_caracf(zcmcf*(i_zone-1)+4) = coef_frot
    v_sdcont_caracf(zcmcf*(i_zone-1)+5) = algo_frot
    v_sdcont_caracf(zcmcf*(i_zone-1)+6) = coef_coul_frot
    v_sdcont_caracf(zcmcf*(i_zone-1)+7) = seuil_init
    v_sdcont_caracf(zcmcf*(i_zone-1)+13) = seuil_auto
    v_sdcont_caracf(zcmcf*(i_zone-1)+16) = glis_maxi
!
! - Contact nodes excluded
!
    if (s_algo_cont .ne. 'LAC') then
        call getvtx(factorKeyword, 'SANS_GROUP_NO', iocc=i_zone, &
                    scal=s_cont_excl, nbret=nb_cont_excl_1)
        call getvtx(factorKeyword, 'SANS_NOEUD', iocc=i_zone, &
                    scal=s_cont_excl, nbret=nb_cont_excl_2)
        l_node_excl = (nb_cont_excl_1 .ne. 0) .or. (nb_cont_excl_2 .ne. 0)
        call getvtx(factorKeyword, 'SANS_GROUP_MA', iocc=i_zone, &
                    scal=s_cont_excl, nbret=nb_cont_excl_3)
        call getvtx(factorKeyword, 'SANS_MAILLE', iocc=i_zone, &
                    scal=s_cont_excl, nbret=nb_cont_excl_4)
        l_node_excl = l_node_excl .or. ((nb_cont_excl_3 .ne. 0) .or. (nb_cont_excl_4 .ne. 0))
    end if
!
! - Friction nodes excluded
!
    if (s_algo_cont .ne. 'LAC') then
        call getvtx(factorKeyword, 'SANS_GROUP_NO_FR', iocc=i_zone, &
                    scal=s_frot_excl, nbret=nb_frot_excl_1)
        call getvtx(factorKeyword, 'SANS_NOEUD_FR', iocc=i_zone, &
                    scal=s_frot_excl, nbret=nb_frot_excl_2)
        l_frot_excl = (nb_frot_excl_1 .ne. 0) .or. (nb_frot_excl_2 .ne. 0)
    end if
!
! - For friction direction to exclude (vector)
!
    if (l_frot_excl) then
        call getvr8(factorKeyword, 'DIRE_EXCL_FROT', iocc=i_zone, nbval=3, vect=dire_excl_frot, &
                    nbret=nb_dire_excl)
        l_dire_excl_frot = (nb_dire_excl .ne. 0)
        if (.not. l_dire_excl_frot) then
! --------- All directions excluded
            dire_excl_frot_i = 2.d0
            dire_excl_frot(1) = 0.d0
            dire_excl_frot(2) = 0.d0
            dire_excl_frot(3) = 0.d0
        else
! --------- Only one direction excluded
            dire_excl_frot_i = 1.d0
        end if
    else
        dire_excl_frot_i = 0.d0
        dire_excl_frot(1) = 0.d0
        dire_excl_frot(2) = 0.d0
        dire_excl_frot(3) = 0.d0
    end if
!
! - Excluded: only node integration scheme
!
    if (.not. l_inte_node) then
        if (l_node_excl .or. l_frot_excl) then
            call utmess('F', 'CONTACT_97')
        end if
        if (.not. mminfl(sdcont_defi, 'MAIT', i_zone)) then
            call utmess('F', 'CONTACT_98')
        end if
    end if

! - Initial contact
    call getvtx(factorKeyword, 'CONTACT_INIT', iocc=i_zone, scal=s_cont_init)
    if (s_cont_init .eq. 'OUI') then
        contInit = 1.d0
    else if (s_cont_init .eq. 'INTERPENETRE') then
        contInit = 2.d0
        contInitDist = -1.d0
        call getvr8(factorKeyword, 'DIST_MAXI', iocc=i_zone, scal=contInitDist, &
                    nbret=nbRet)
        if (nbRet == 0) then
            contInitDist = -1.d0
        end if
    else if (s_cont_init .eq. 'NON') then
        contInit = 0.d0
    else
        ASSERT(ASTER_FALSE)
    end if

!
! - Bilateral contact
!
    if (s_algo_cont .ne. 'LAC') then
        call getvtx(factorKeyword, 'GLISSIERE', iocc=i_zone, scal=s_gliss)
        if (s_gliss .eq. 'OUI') then
            l_gliss = ASTER_TRUE
        else if (s_gliss .eq. 'NON') then
            l_gliss = ASTER_FALSE
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
! - Large sliding method
!
    if (l_frot .and. l_cont_cont) then
        call getvtx(factorKeyword, 'GRAND_GLIS', iocc=i_zone, scal=s_grglis, nbret=noc)
        if (noc .ge. 1) then
            if (s_grglis(1:3) .eq. 'OUI') then
                l_granglis = ASTER_TRUE
            elseif (s_grglis(1:3) .eq. 'NON') then
                l_granglis = ASTER_FALSE
            else
                print s_grglis(1:3)
                ASSERT(ASTER_FALSE)
            end if
        else
            l_granglis = ASTER_FALSE
        end if
        if (l_granglis) then
            v_sdcont_caracf(zcmcf*(i_zone-1)+15) = 1.d0
        else
            v_sdcont_caracf(zcmcf*(i_zone-1)+15) = 0.d0
        end if
    end if
!
    v_sdcont_caracf(zcmcf*(i_zone-1)+8) = contInit
    v_sdcont_caracf(zcmcf*(i_zone-1)+17) = contInitDist
    if (l_gliss) then
        v_sdcont_caracf(zcmcf*(i_zone-1)+9) = 1.d0
    else
        v_sdcont_caracf(zcmcf*(i_zone-1)+9) = 0.d0
    end if
    if (l_node_excl) then
        v_sdcont_caracf(zcmcf*(i_zone-1)+10) = 1.d0
    else
        v_sdcont_caracf(zcmcf*(i_zone-1)+10) = 0.d0
    end if
    if (l_frot_excl) then
        v_sdcont_caracf(zcmcf*(i_zone-1)+11) = 1.d0
    else
        v_sdcont_caracf(zcmcf*(i_zone-1)+11) = 0.d0
    end if
    v_sdcont_caracf(zcmcf*(i_zone-1)+12) = dire_excl_frot_i
    v_sdcont_exclfr(zexcl*(i_zone-1)+1) = dire_excl_frot(1)
    v_sdcont_exclfr(zexcl*(i_zone-1)+2) = dire_excl_frot(2)
    v_sdcont_exclfr(zexcl*(i_zone-1)+3) = dire_excl_frot(3)
!
end subroutine
