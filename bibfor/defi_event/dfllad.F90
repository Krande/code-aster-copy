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
subroutine dfllad(sdlist)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dinogd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/reliem.h"
#include "asterfort/utcmp2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "event_def.h"
!
    character(len=8), intent(in) :: sdlist
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_LIST_INST - Read parameters
!
! Keyword ADAPTATION
!
! --------------------------------------------------------------------------------------------------
!
! In  sdlist           : name of DEFI_LIST_INST datastructure
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'ADAPTATION'
    integer(kind=8) :: nbAdap, nbret
    integer(kind=8) :: ibid, nb_iter_newton_ref, nb_incr_seuil
    integer(kind=8) :: iAdap
    character(len=16) :: event_typek, nom_para, crit_comp, action_typek, nom_cham
    character(len=8) :: nomgd, nom_cmp
    real(kind=8) :: pcent_augm, vale_ref, valer
    integer(kind=8) :: valei, nucmp(1)
    character(len=24) :: sdlistAEvenrName
    real(kind=8), pointer :: sdlistAEvenr(:) => null()
    character(len=24) :: sdlistALocaName
    integer(kind=8), pointer :: sdlistALoca(:) => null()
    integer(kind=8), pointer :: v_lst_loca(:) => null()
    character(len=24) :: sdlistATplurName
    real(kind=8), pointer :: sdlistATplur(:) => null()
    character(len=24) :: sdlistATplukName
    character(len=16), pointer :: sdlistATpluk(:) => null()
    character(len=24) :: sdlist_linfor
    real(kind=8), pointer :: v_sdlist_linfor(:) => null()
    integer(kind=8) :: nocc, nb_loca, lg_ini, etat_loca
    character(len=8)  :: mesh
    character(len=24) :: model, lst_loca
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!

    nbAdap = 0
!
! - Access to datastructures
!
    sdlist_linfor = sdlist(1:8)//'.LIST.INFOR'
    call jeveuo(sdlist_linfor, 'E', vr=v_sdlist_linfor)

! - Get number of adaptation keywords
    call getfac(factorKeyword, nbAdap)

! - For IMPLEX: only one MODE_CALCUL_TPLUS
    call getvtx('ADAPTATION', 'MODE_CALCUL_TPLUS', iocc=nbAdap, scal=action_typek, &
                nbret=nbret)
    if (nbAdap .ne. 1 .and. action_typek .eq. adapActionKeyword(ADAP_ACT_IMPLEX)) then
        call utmess('F', 'DISCRETISATION_15')
    end if

! - Create datastructure
    sdlistAEvenrName = sdlist(1:8)//'.ADAP.EVENR'
    sdlistALocaName = sdlist(1:8)//'.ADAP.LOCA'
    sdlistATplurName = sdlist(1:8)//'.ADAP.TPLUR'
    sdlistATplukName = sdlist(1:8)//'.ADAP.TPLUK'
    call wkvect(sdlistAEvenrName, 'G V R', nbAdap*SIZE_LAEVR, vr=sdlistAEvenr)
    call wkvect(sdlistALocaName, 'G V I', nbAdap*SIZE_LALOCA, vi=sdlistALoca)
    call wkvect(sdlistATplurName, 'G V R', nbAdap*SIZE_LATPR, vr=sdlistATplur)
    call wkvect(sdlistATplukName, 'G V K16', nbAdap*SIZE_LATPK, vk16=sdlistATpluk)
    v_sdlist_linfor(10) = nbAdap
    sdlistALoca(1:nbAdap*SIZE_LALOCA) = 0

! - Read parameters
    do iAdap = 1, nbAdap
! ----- Get event
        call getvtx(factorKeyword, 'EVENEMENT', iocc=iAdap, scal=event_typek, nbret=nbret)
        if (event_typek .eq. adapEventKeyword(ADAP_EVT_NONE)) then
            sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+1) = ADAP_EVT_NONE
            call utmess('A', 'DISCRETISATION_5')
        else if (event_typek .eq. adapEventKeyword(ADAP_EVT_ALLSTEPS)) then
            sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+1) = ADAP_EVT_ALLSTEPS
        else if (event_typek .eq. adapEventKeyword(ADAP_EVT_TRIGGER)) then
            sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+1) = ADAP_EVT_TRIGGER
        else
            ASSERT(.false.)
        end if
!
! ----- Options for 'SEUIL'
!
        if (event_typek .eq. adapEventKeyword(ADAP_EVT_TRIGGER)) then
            call getvis(factorKeyword, 'NB_INCR_SEUIL', iocc=iAdap, scal=nb_incr_seuil, nbret=nbret)
            sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+2) = nb_incr_seuil
            call getvtx(factorKeyword, 'NOM_PARA', iocc=iAdap, scal=nom_para, nbret=nbret)
            if (nom_para .eq. 'NB_ITER_NEWTON') then
                sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+3) = 1.d0
                valei = 0
                call getvis(factorKeyword, 'VALE_I', iocc=iAdap, scal=valei, nbret=nbret)
                valer = valei
            else
                ASSERT(.false.)
            end if
            sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+5) = valer
            call getvtx(factorKeyword, 'CRIT_COMP', iocc=iAdap, scal=crit_comp, nbret=nbret)
            if (crit_comp .eq. 'LT') then
                sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+4) = 1.d0
            else if (crit_comp .eq. 'GT') then
                sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+4) = 2.d0
            else if (crit_comp .eq. 'LE') then
                sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+4) = 3.d0
            else if (crit_comp .eq. 'GE') then
                sdlistAEvenr(SIZE_LAEVR*(iAdap-1)+4) = 4.d0
            else
                ASSERT(.false.)
            end if
        end if
!
! ----- Options for 'MODE_CALCUL_TPLUS'
!
        call getvtx(factorKeyword, 'MODE_CALCUL_TPLUS', iocc=iAdap, scal=action_typek, nbret=nbret)
        if (action_typek .eq. adapActionKeyword(ADAP_ACT_FIXE)) then
            sdlistATplur(SIZE_LATPR*(iAdap-1)+1) = ADAP_ACT_FIXE
        else if (action_typek .eq. adapActionKeyword(ADAP_ACT_INCR_QUANT)) then
            sdlistATplur(SIZE_LATPR*(iAdap-1)+1) = ADAP_ACT_INCR_QUANT
        else if (action_typek .eq. adapActionKeyword(ADAP_ACT_ITER)) then
            sdlistATplur(SIZE_LATPR*(iAdap-1)+1) = ADAP_ACT_ITER
        else if (action_typek .eq. adapActionKeyword(ADAP_ACT_IMPLEX)) then
            sdlistATplur(SIZE_LATPR*(iAdap-1)+1) = ADAP_ACT_IMPLEX
            if (event_typek .ne. adapEventKeyword(ADAP_EVT_ALLSTEPS)) then
                call utmess('F', 'DISCRETISATION_14')
            end if
        else
            ASSERT(.false.)
        end if
!
! ----- Options for MODE_CALCUL_TPLUS/FIXE
!
        if (action_typek .eq. adapActionKeyword(ADAP_ACT_FIXE)) then
            call getvr8(factorKeyword, 'PCENT_AUGM', iocc=iAdap, scal=pcent_augm, nbret=nbret)
            sdlistATplur(SIZE_LATPR*(iAdap-1)+2) = pcent_augm
        end if
!
! ----- Options for MODE_CALCUL_TPLUS/DELTA_GRANDEUR
!
        if (action_typek .eq. adapActionKeyword(ADAP_ACT_INCR_QUANT)) then
            call getvr8(factorKeyword, 'VALE_REF', iocc=iAdap, scal=vale_ref, nbret=nbret)
            sdlistATplur(SIZE_LATPR*(iAdap-1)+3) = vale_ref
            call getvtx(factorKeyword, 'NOM_PARA', iocc=iAdap, scal=nom_para, nbret=nbret)
            call getvtx(factorKeyword, 'NOM_CHAM', iocc=iAdap, scal=nom_cham, nbret=nbret)
            call getvtx(factorKeyword, 'NOM_CMP', iocc=iAdap, scal=nom_cmp, nbret=nbret)

            if (nom_cham .eq. 'DEPL') then
                call getvtx(factorKeyword, 'GROUP_NO', iocc=iAdap, nbret=nocc)
            else if (nom_cham .eq. 'SIEF_ELGA' .or. nom_cham .eq. 'VARI_ELGA') then
                call getvtx(factorKeyword, 'GROUP_MA', iocc=iAdap, nbret=nocc)
            else
                ASSERT(.false.)
            end if
            etat_loca = merge(LOCA_TOUT, LOCA_PARTIEL, nocc .eq. 0)

            if (etat_loca .eq. LOCA_PARTIEL) then
                call getvid(' ', 'MODELE', scal=model, nbret=nocc)
                if (nocc .ne. 1) call utmess('F', 'LISTINST_4')
                call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

                write (lst_loca, '(A19,I5.5)') '&&OP0028.ADAP.LOCA.', iAdap
                if (nom_cham .eq. 'DEPL') then
                    call reliem(model, mesh, 'NU_NOEUD', factorKeyword, iAdap, 1, ['GROUP_NO'], &
                                ['GROUP_NO'], lst_loca, nb_loca)
                else if (nom_cham .eq. 'SIEF_ELGA' .or. nom_cham .eq. 'VARI_ELGA') then
                    call reliem(model, mesh, 'NU_MAILLE', factorKeyword, iAdap, 1, ['GROUP_MA'], &
                                ['GROUP_MA'], lst_loca, nb_loca)
                else
                    ASSERT(.false.)
                end if
                if (nb_loca .eq. 0) etat_loca = LOCA_VIDE
            end if

            if (etat_loca .eq. LOCA_VIDE) then
                call jeveuo(sdlistALocaName, 'E', vi=sdlistALoca)
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+1) = etat_loca
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+2) = 0
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+3) = 0
            else if (etat_loca .eq. LOCA_PARTIEL) then
                call jelira(lst_loca, 'LONMAX', nb_loca)
                call jelira(sdlistALocaName, 'LONMAX', lg_ini)
                call juveca(sdlistALocaName, lg_ini+nb_loca)

                call jeveuo(sdlistALocaName, 'E', vi=sdlistALoca)
                call jeveuo(lst_loca, 'L', vi=v_lst_loca)

                sdlistALoca(SIZE_LALOCA*(iAdap-1)+1) = etat_loca
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+2) = lg_ini+1
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+3) = lg_ini+nb_loca
                sdlistALoca(lg_ini+1:lg_ini+nb_loca) = v_lst_loca(1:nb_loca)
            else if (etat_loca .eq. LOCA_TOUT) then
                call jeveuo(sdlistALocaName, 'E', vi=sdlistALoca)
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+1) = etat_loca
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+2) = 0
                sdlistALoca(SIZE_LALOCA*(iAdap-1)+3) = 0
            end if

            nomgd = dinogd(nom_cham)
            call utcmp2(nomgd, factorKeyword, iAdap, 1, nom_cmp, &
                        nucmp, ibid)
            sdlistATpluk(SIZE_LATPK*(iAdap-1)+1) = nom_para
            sdlistATpluk(SIZE_LATPK*(iAdap-1)+2) = nom_cham
            sdlistATpluk(SIZE_LATPK*(iAdap-1)+3) = nom_cmp
            sdlistATplur(SIZE_LATPR*(iAdap-1)+4) = nucmp(1)
        end if
!
! ----- Options for MODE_CALCUL_TPLUS/ITER_NEWTON
!
        if (action_typek .eq. adapActionKeyword(ADAP_ACT_ITER)) then
            call getvis(factorKeyword, 'NB_ITER_NEWTON_REF', iocc=iAdap, scal=nb_iter_newton_ref, &
                        nbret=nbret)
            sdlistATplur(SIZE_LATPR*(iAdap-1)+5) = nb_iter_newton_ref
        end if
    end do
!
    call jedema()
end subroutine
