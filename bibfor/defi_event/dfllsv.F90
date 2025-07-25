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
! aslint: disable=W1504
!
subroutine dfllsv(v_sdlist_linfor, v_sdlist_eevenr, v_sdlist_eevenk, sdlist_loca, &
                  v_sdlist_esubdr, i_fail_save, &
                  event_typek, vale_ref, nom_cham, nom_cmp, &
                  crit_cmp, lst_loca, etat_loca, pene_maxi, resi_glob_maxi, &
                  action_typek, subd_methode, subd_auto, subd_pas_mini, &
                  subd_pas, subd_niveau, pcent_iter_plus, coef_maxi, &
                  subd_inst, subd_duree)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
!
    real(kind=8), pointer :: v_sdlist_linfor(:)
    real(kind=8), pointer :: v_sdlist_eevenr(:)
    character(len=16), pointer :: v_sdlist_eevenk(:)
    character(len=24), intent(in):: sdlist_loca
    real(kind=8), pointer :: v_sdlist_esubdr(:)
    integer(kind=8), intent(in) :: i_fail_save
    character(len=16), intent(in) :: event_typek
    real(kind=8), intent(in) :: vale_ref
    character(len=16), intent(in) :: nom_cham
    character(len=16), intent(in) :: nom_cmp
    character(len=16), intent(in) :: crit_cmp
    character(len=24), intent(in) :: lst_loca
    integer(kind=8), intent(in):: etat_loca
    real(kind=8), intent(in) :: pene_maxi
    real(kind=8), intent(in) :: resi_glob_maxi
    character(len=16), intent(in) :: action_typek
    character(len=16), intent(in) :: subd_methode
    real(kind=8), intent(in) :: subd_pas_mini
    integer(kind=8), intent(in) :: subd_niveau
    integer(kind=8), intent(in) :: subd_pas
    character(len=16), intent(in) :: subd_auto
    real(kind=8), intent(in) :: subd_inst
    real(kind=8), intent(in) :: subd_duree
    real(kind=8), intent(in) :: pcent_iter_plus
    real(kind=8), intent(in) :: coef_maxi
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_LIST_INST - Read parameters
!
! Save parameters in datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  v_sdlist_linfor  : pointer to LIST.INFOR object
! In  v_sdlist_eevenr  : pointer to ECHEC.EVENR object
! In  v_sdlist_eevenk  : pointer to ECHEC.EVENK object
! In  v_sdlist_esubdr  : pointer to ECHEC.SUBDR object
! In  i_fail_save      : current index for ECHEC keyword
! In  event_typek      : type of event
! In  vale_ref         : value of VALE_REF for EVENEMENT=DELTA_GRANDEUR
! In  nom_cham         : value of NOM_CHAM for EVENEMENT=DELTA_GRANDEUR
! In  nom_cmp          : value of NOM_CMP for EVENEMENT=DELTA_GRANDEUR
! In  crit_cmp         : value of CRIT_CMP for EVENEMENT=DELTA_GRANDEUR
! In  lst_loca         : list of cells or nodes for EVENEMENT=DELTA_GRANDEUR
! In  etat_loca        : in relation with lst_loca 0=void, 1=partial (lst_loca), 2=all
! In  pene_maxi        : value of PENE_MAXI for EVENEMENT=INTERPENETRATION
! In  resi_glob_maxi   : value of RESI_GLOB_MAXI for EVENEMENT=RESI_MAXI
! In  action_typek     : type of action
! In  subd_methode     : value of SUBD_METHODE for ACTION=DECOUPE
! In  subd_pas_mini    : value of SUBD_PAS_MINI for ACTION=DECOUPE
! In  subd_niveau      : value of SUBD_NIVEAU for ACTION=DECOUPE
! In  subd_pas         : value of SUBD_PAS for ACTION=DECOUPE
! In  subd_auto        : value of SUBD_AUTO for ACTION=DECOUPE
! In  subd_inst        : value of SUBD_INST for ACTION=DECOUPE
! In  subd_duree       : value of SUBD_DUREE for ACTION=DECOUPE
! In  pcent_iter_plus  : value of PCENT_ITER_PLUS for ACTION=ITER_SUPPL
! In  coef_maxi        : value of COEF_MAXI for ACTION=ADAPT_COEF_PENA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8):: nb_loca, lg_ini
    integer(kind=8), pointer :: v_sdlist_loca(:) => null()
    integer(kind=8), pointer :: v_lst_loca(:) => null()
! --------------------------------------------------------------------------------------------------

    call jemarq()
!
! - Alarm for no-step cut
!
    if (event_typek .eq. failActionKeyword(FAIL_EVT_ERROR)) then
        if (action_typek .eq. failEventKeyword(FAIL_ACT_STOP)) then
            call utmess('I', 'DISCRETISATION_9')
        end if
    end if
!
! - At least one ACTION=DECOUPE
!
    if (action_typek .eq. failActionKeyword(FAIL_ACT_CUT)) then
        v_sdlist_linfor(7) = 1.d0
    end if
!
! - Type of event
!
    if (event_typek .eq. failEventKeyword(FAIL_EVT_ERROR)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_ERROR
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_INCR_QUANT)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_INCR_QUANT
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_COLLISION)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_COLLISION
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_INTERPENE)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_INTERPENE
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_DIVE_RESI)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_DIVE_RESI
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_INSTABILITY)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_INSTABILITY
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_RESI_MAXI)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+1) = FAIL_EVT_RESI_MAXI
    else
        ASSERT(.false.)
    end if
!
! - Type of action
!
    if (action_typek .eq. failActionKeyword(FAIL_ACT_STOP)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+2) = FAIL_ACT_STOP
    else if (action_typek .eq. failActionKeyword(FAIL_ACT_CUT)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+2) = FAIL_ACT_CUT
    else if (action_typek .eq. failActionKeyword(FAIL_ACT_ITER)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+2) = FAIL_ACT_ITER
    else if (action_typek .eq. failActionKeyword(FAIL_ACT_ADAPT_COEF)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+2) = FAIL_ACT_ADAPT_COEF
    else if (action_typek .eq. failActionKeyword(FAIL_ACT_CONTINUE)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+2) = FAIL_ACT_CONTINUE
    else
        ASSERT(.false.)
    end if
!
! - Parameters for EVENEMENT = 'DELTA_GRANDEUR'
!
    if (event_typek .eq. failEventKeyword(FAIL_EVT_INCR_QUANT)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+5) = vale_ref
        v_sdlist_eevenk(SIZE_LEEVK*(i_fail_save-1)+1) = nom_cham
        v_sdlist_eevenk(SIZE_LEEVK*(i_fail_save-1)+2) = nom_cmp
        v_sdlist_eevenk(SIZE_LEEVK*(i_fail_save-1)+3) = crit_cmp

        if (etat_loca .eq. LOCA_VIDE) then
            call jeveuo(sdlist_loca, 'E', vi=v_sdlist_loca)
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+1) = etat_loca
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+2) = 0
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+3) = 0
        else if (etat_loca .eq. LOCA_PARTIEL) then
            call jelira(lst_loca, 'LONMAX', nb_loca)
            call jelira(sdlist_loca, 'LONMAX', lg_ini)
            call juveca(sdlist_loca, lg_ini+nb_loca)

            call jeveuo(sdlist_loca, 'E', vi=v_sdlist_loca)
            call jeveuo(lst_loca, 'L', vi=v_lst_loca)

            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+1) = etat_loca
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+2) = lg_ini+1
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+3) = lg_ini+nb_loca
            v_sdlist_loca(lg_ini+1:lg_ini+nb_loca) = v_lst_loca(1:nb_loca)
        else if (etat_loca .eq. LOCA_TOUT) then
            call jeveuo(sdlist_loca, 'E', vi=v_sdlist_loca)
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+1) = etat_loca
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+2) = 0
            v_sdlist_loca(SIZE_LELOCA*(i_fail_save-1)+3) = 0
        end if
    end if
!
! - Parameters for EVENEMENT = 'INTERPENETRATION'
!
    if (event_typek .eq. failEventKeyword(FAIL_EVT_INTERPENE)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+6) = pene_maxi
    end if
!
! - Parameters for EVENEMENT = 'RESI_MAXI'
!
    if (event_typek .eq. failEventKeyword(FAIL_EVT_RESI_MAXI)) then
        v_sdlist_eevenr(SIZE_LEEVR*(i_fail_save-1)+7) = resi_glob_maxi
    end if
!
! - Parameters for ACTION = 'DECOUPE'
!
    if (action_typek .ne. failActionKeyword(FAIL_ACT_STOP)) then
        if (subd_methode .eq. 'MANUEL') then
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+1) = 1.d0
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+2) = subd_pas
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+3) = subd_pas_mini
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+4) = subd_niveau
        else if (subd_methode .eq. 'AUTO') then
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+1) = 2.d0
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+3) = subd_pas_mini
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+5) = subd_inst
            v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+6) = subd_duree
            if (subd_auto .eq. 'COLLISION') then
                v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+10) = 1.d0
            else if (subd_auto .eq. 'EXTRAPOLE') then
                v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+10) = 2.d0
            else
                ASSERT(.false.)
            end if
        end if
    end if
!
! - Parameters for ACTION = 'ITER_SUPPL'
!
    if (action_typek .eq. failActionKeyword(FAIL_ACT_ITER)) then
        v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+7) = pcent_iter_plus
    end if
!
! - Parameters for ACTION = 'ADAPT_COEF_PENA'
!
    if (action_typek .eq. failActionKeyword(FAIL_ACT_ADAPT_COEF)) then
        v_sdlist_esubdr(SIZE_LESUR*(i_fail_save-1)+8) = coef_maxi
    end if
!
    call jedema()
!
end subroutine
