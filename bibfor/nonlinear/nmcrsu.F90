! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine nmcrsu(sddisc, listInst, ds_conv, &
                  l_implex_, ds_contact_)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getAdapAction.h"
#include "asterfort/getAdapEvent.h"
#include "asterfort/GetResi.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedup1.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcerr.h"
#include "asterfort/nmcrld.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "event_def.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc, listInst
    type(NL_DS_Conv), intent(in) :: ds_conv
    aster_logical, optional, intent(in) :: l_implex_
    type(NL_DS_Contact), optional, intent(in) :: ds_contact_
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (STRUCTURES DE DONNES)
!
! CREATION SD DISCRETISATION - SUBDIVISION AUTO
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_conv          : datastructure for convergence management
! In  ds_algopara      : datastructure for algorithm parameters
! IN  LISINS : SD_LIST_INST OU SD_LISTR8
! In  sddisc           : datastructure for time discretization
! IN  LIMPEX : .TRUE. SI IMPLEX
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: metlis
    integer(kind=8) :: iret
    real(kind=8) :: pas_mini_elas, valr
    integer(kind=8) :: nbAdap, iAdap
    integer(kind=8) :: iter_glob_maxi, iter_glob_elas
    integer(kind=8) :: ifm, niv, itmx, vali
    aster_logical :: l_implex
    real(kind=8) :: resi_glob_maxi, resi_glob_rela, inikry
    character(len=16) :: dsType, nopara
    character(len=24) :: lisevr, lisevk, liseloca, lisesu
    character(len=24) :: lisavr, lisaloca, listpr, listpk
    character(len=24) :: tpsevr, tpsevk, tpseloca, tpsesu
    character(len=24) :: tpsavr, tpsaloca, tpstpr, tpstpk
    character(len=24) :: tpsext
    integer(kind=8) :: jtpsex, eventType, action_type
    character(len=8), pointer:: v_modele_dli(:) => null()
    character(len=8):: modele_dli, modele_snl
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_16')
    end if

! - Get parameters for convergence
    call GetResi(ds_conv, type='RESI_GLOB_RELA', user_para_=resi_glob_rela)
    call GetResi(ds_conv, type='RESI_GLOB_MAXI', user_para_=resi_glob_maxi)
    iter_glob_maxi = ds_conv%iter_glob_maxi
    iter_glob_elas = ds_conv%iter_glob_elas

! - Initializations
    l_implex = ASTER_FALSE
    if (present(l_implex_)) then
        l_implex = l_implex_
    end if
    inikry = 0.9d0
    pas_mini_elas = 0.d0
    call utdidt('L', sddisc, 'LIST', 'NADAPT', vali_=nbAdap)
    call utdidt('L', sddisc, 'LIST', 'METHODE', valk_=metlis)
!
! --- NOM SDS DE LA LISINS
!
    lisevr = listInst(1:8)//'.ECHE.EVENR'
    lisevk = listInst(1:8)//'.ECHE.EVENK'
    liseloca = listInst(1:8)//'.ECHE.LOCA'
    lisesu = listInst(1:8)//'.ECHE.SUBDR'
    lisavr = listInst(1:8)//'.ADAP.EVENR'
    lisaloca = listInst(1:8)//'.ADAP.LOCA'
    listpr = listInst(1:8)//'.ADAP.TPLUR'
    listpk = listInst(1:8)//'.ADAP.TPLUK'
!
! --- NOM SDS DE LA SDDISC
!
    tpsevr = sddisc(1:19)//'.EEVR'
    tpsevk = sddisc(1:19)//'.EEVK'
    tpseloca = sddisc(1:19)//'.ELOC'
    tpsesu = sddisc(1:19)//'.ESUR'
    tpsavr = sddisc(1:19)//'.AEVR'
    tpsaloca = sddisc(1:19)//'.ALOC'
    tpstpr = sddisc(1:19)//'.ATPR'
    tpstpk = sddisc(1:19)//'.ATPK'
!
! --- LECTURE DE LA LISTE D'INSTANTS
!
    call gettco(listInst, dsType)
!
    if (dsType .eq. 'LISTR8_SDASTER') then
!
! ----- CREATION EVENEMENTS ERREURS: ARRET
!
        call nmcrld(sddisc)
    else if (dsType .eq. 'LIST_INST') then
!
! ----- COPIE LOCALE DES OBJETS DE LA LISINS
!
        call jedup1(lisevr, 'V', tpsevr)
        call jedup1(lisevk, 'V', tpsevk)
        call jedup1(liseloca, 'V', tpseloca)
        call jedup1(lisesu, 'V', tpsesu)
        if (nbAdap .ne. 0) then
            call jedup1(lisavr, 'V', tpsavr)
            call jedup1(lisaloca, 'V', tpsaloca)
            call jedup1(listpr, 'V', tpstpr)
            call jedup1(listpk, 'V', tpstpk)
        end if

        ! Verification de la coherence du modele si renseigne dans DEFI_LIST_INST
        call jeveuo(listInst(1:8)//'.MODELE', 'L', vk8=v_modele_dli)
        modele_dli = v_modele_dli(1)
        if (modele_dli .ne. ' ') then
            call getvid(' ', 'MODELE', scal=modele_snl)
            if (modele_dli .ne. modele_snl) call utmess('F', 'MECANONLINE_7')

        end if
    else
        ASSERT(ASTER_FALSE)
    end if

!
! --- EN GESTION AUTO, AVEC UN CRITERE D'ADAPTATION EN SEUIL SUR
!     NB_ITER_NEWT, ON MET VALE = ITER_GLOB_MAXI/2 SI VALE N'A PAS
!     ETE RENSIGNE DANS DEFI_LIST_INST
!     ON NE CONSIDERE PAS LE CAS DE DE ITER_GLOB_ELAS CAR C'EST ACTIVE
!     (MATRICE SECANTE) QU'EN CAS DE DIFFICULTE
!
    if (metlis .eq. 'AUTO') then
        call getvis('CONVERGENCE', 'ITER_GLOB_MAXI', iocc=1, scal=itmx, nbret=iret)
        do iAdap = 1, nbAdap
            call getAdapEvent(sddisc, iAdap, eventType)
            if (eventType .eq. ADAP_EVT_TRIGGER) then
                call utdidt('L', sddisc, 'ADAP', 'NOM_PARA', index_=iAdap, &
                            valk_=nopara)
                if (nopara .eq. 'NB_ITER_NEWT') then
                    call utdidt('L', sddisc, 'ADAP', 'VALE', index_=iAdap, &
                                vali_=vali)
                    if (vali .eq. 0) then
                        vali = itmx/2
                        valr = vali
                        call utdidt('E', sddisc, 'ADAP', 'VALE', index_=iAdap, &
                                    valr_=valr)
                    end if
                end if
            end if
        end do
    end if

! - VERIF COHERENCE AVEC IMPLEX
    if (metlis .eq. 'AUTO') then
        do iAdap = 1, nbAdap
            call getAdapAction(sddisc, iAdap, action_type)
            if (action_type .eq. ADAP_ACT_IMPLEX) then
                if (.not. l_implex) then
                    call utmess('F', 'MECANONLINE6_4')
                end if
            end if
        end do
    end if

! - CREATION SD STOCKAGE DES INFOS EN COURS DE CALCUL
    if (present(ds_contact_)) then
        call nmcerr(sddisc, iter_glob_maxi, iter_glob_elas, pas_mini_elas, resi_glob_maxi, &
                    resi_glob_rela, inikry, ds_contact_)
    else
        call nmcerr(sddisc, iter_glob_maxi, iter_glob_elas, pas_mini_elas, resi_glob_maxi, &
                    resi_glob_rela, inikry)
    end if
!
! --- OBJET POUR PROLONGEMENT DECOUPE
!
    tpsext = sddisc(1:19)//'.AEXT'
    call wkvect(tpsext, 'V V R', 3, jtpsex)
    zr(jtpsex-1+1) = r8vide()
    zr(jtpsex-1+2) = r8vide()
    zr(jtpsex-1+3) = r8vide()
!
    call jedema()
!
end subroutine
