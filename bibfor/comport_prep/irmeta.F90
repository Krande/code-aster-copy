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
subroutine irmeta(ifi, field_med, meta_elno, field_loca, ligrel, &
                  nb_cmp_sele, cmp_name_sele, partie, numpt, instan, &
                  nume_store, nbmaec, limaec, result, cara_elem, &
                  lfichUniq, codret)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescrm.h"
#include "asterfort/cesexi.h"
#include "asterfort/codent.h"
#include "asterfort/comp_meca_uvar.h"
#include "asterfort/comp_meta_pvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/irceme.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8), intent(in) :: ifi
    character(len=64), intent(in) :: field_med
    character(len=19), intent(in) :: meta_elno
    character(len=8), intent(in) :: field_loca
    character(len=19), intent(in) :: ligrel
    integer(kind=8), intent(in) :: nb_cmp_sele
    character(len=*), intent(in) :: cmp_name_sele(*)
    character(len=*), intent(in) :: partie
    integer(kind=8), intent(in) :: numpt
    real(kind=8), intent(in) :: instan
    integer(kind=8), intent(in) :: nume_store
    integer(kind=8), intent(in) :: nbmaec
    integer(kind=8), intent(in) :: limaec(*)
    character(len=8), intent(in) :: result
    character(len=8), intent(in) :: cara_elem
    aster_logical, intent(in) :: lfichUniq
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Post-treatment (IMPR_RESU)
!
! Create META_ELNO_NOMME for name of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_store       : index to store in results
! In  field_med        : name of MED field
! In  field_loca       : localization of field
!                        /'ELNO'/'ELGA'/'ELEM'
! In  result           : name of results datastructure
! In  ligrel           : name of ligrel of the field
! Out codret           : error code
!                        0   - Everything is OK
!                        400 - error
!       IFI    : UNITE LOGIQUE D'IMPRESSION DU CHAMP
!       PARTIE : IMPRESSION DE LA PARTIE IMAGINAIRE OU REELLE POUR
!                UN CHAMP COMPLEXE
!       nb_cmp_sele  : NOMBRE DE COMPOSANTES A ECRIRE
!       cmp_name_sele : NOMS DES COMPOSANTES A ECRIRE
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       INSTAN : VALEUR DE L'INSTANT A ARCHIVER
!       NUMORD : NUMERO D'ORDRE DU CHAMP
!       NBMAEC : NOMBRE DE MAILLES A ECRIRE (0, SI TOUTES LES MAILLES)
!       LIMAEC : LISTE DES MAILLES A ECRIRE SI EXTRAIT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: mapZoneNume, i_elem, i_pt, i_vari, i_vari_redu, i_spt
    integer(kind=8) :: nb_vari, nb_pt, nb_spt, nb_vari_zone
    integer(kind=8) :: nb_vari_redu, mapNbZone, nb_elem, nbVariMaxi, nbCellMesh, nbCellInZone
    integer(kind=8) :: ntVari, codret_dummy
    integer(kind=8) :: posit, iret, affe_type, affe_indx, nume_elem, nbCmpDyna
    integer(kind=8) :: jv_elno_cesd, jv_elno_cesl, jv_elnr_cesd, jv_elnr_cesl, jv_elno, jv_elnr
    character(len=7) :: saux07
    character(len=8) :: saux08
    character(len=8), parameter :: base_name = '&&IRMETA'
    character(len=19) :: comporMeta
    character(len=19), parameter :: meta_elno_s = '&&IRMETA.VARIELGA_S'
    character(len=19), parameter :: meta_elnr = '&&IRMETA.VARIELGR'
    character(len=19), parameter :: meta_elnr_s = '&&IRMETA.VARIELGR_S'
    character(len=19) :: vari_link
    character(len=19), parameter :: vari_redu = '&&IRMETA.VARIREDU'
    integer(kind=8), pointer :: v_vari_link(:) => null()
    character(len=16), pointer :: v_vari_redu(:) => null()
    character(len=19), parameter :: label_med = '&&IRMETA.LABELMED'
    character(len=19), parameter :: label_vxx = '&&IRMETA.LABELVXX'
    character(len=8), pointer :: v_label_vxx(:) => null()
    character(len=16), pointer :: v_label_med(:) => null()
    character(len=64) :: nomres
    real(kind=8), pointer :: v_elnr_cesv(:) => null()
    real(kind=8), pointer :: v_elga_cesv(:) => null()
    integer(kind=8), pointer :: v_compor_desc(:) => null()
    integer(kind=8), pointer :: v_compor_lima(:) => null()
    integer(kind=8), pointer :: v_compor_lima_lc(:) => null()
    character(len=19), parameter :: comporInfo = '&&IRMETA.INFO'
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    character(len=16) :: field_type
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    field_type = 'META_ELNO'
    ASSERT(field_loca .eq. 'ELNO')
    codret = 0

! - Get name of <CARTE> COMPORMETA
    call rsexch(' ', result, 'COMPORMETA', nume_store, comporMeta, iret)
    if (iret .ne. 0) then
        codret = 400
        goto 99
    end if

! - Prepare informations about internal variables
    call comp_meta_pvar(ligrel, comporMeta, comporInfo)

! - Access to informations
    call jeveuo(comporInfo(1:19)//'.INFO', 'L', vi=comporInfoInfo)
    call jeveuo(comporInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)
    nbCellMesh = comporInfoInfo(1)
    mapNbZone = comporInfoInfo(2)
    nbVariMaxi = comporInfoInfo(3)
    ntVari = comporInfoInfo(4)

! - Create list of internal variables and link to zone in <CARTE> COMPOR
    call comp_meca_uvar(comporInfo, base_name, vari_redu, nb_vari_redu, codret)
    call jeveuo(vari_redu, 'L', vk16=v_vari_redu)

! - Access to <CARTE> COMPOR
    call jeveuo(comporMeta//'.DESC', 'L', vi=v_compor_desc)
    call jeveuo(jexnum(comporMeta//'.LIMA', 1), 'L', vi=v_compor_lima)
    call jeveuo(jexatr(comporMeta//'.LIMA', 'LONCUM'), 'L', vi=v_compor_lima_lc)

! - Transform META_ELNO in META_ELNO_S
    call celces(meta_elno, 'V', meta_elno_s)
    call jeveuo(meta_elno_s//'.CESD', 'L', jv_elno_cesd)
    call jeveuo(meta_elno_s//'.CESL', 'L', jv_elno_cesl)
    call jeveuo(meta_elno_s//'.CESV', 'L', vr=v_elga_cesv)
!
! - Prepare objects to reduced list of internal variables
!
    call wkvect(label_vxx, 'V V K8', nb_vari_redu, vk8=v_label_vxx)
    call wkvect(label_med, 'V V K16', 2*nb_vari_redu, vk16=v_label_med)
    do i_vari_redu = 1, nb_vari_redu
        call codent(i_vari_redu, 'G', saux07)
        v_label_vxx(i_vari_redu) = 'V'//saux07
        v_label_med(2*(i_vari_redu-1)+1) = 'V'//saux07
        v_label_med(2*(i_vari_redu-1)+2) = v_vari_redu(i_vari_redu)
    end do
!
! - Create META_ELNR_S on reduced list of internal variables
!
    call cescrm('V', meta_elnr_s, field_loca, 'VARI_R', nb_vari_redu, &
                v_label_vxx, meta_elno_s)
    call jeveuo(meta_elnr_s//'.CESD', 'L', jv_elnr_cesd)
    call jeveuo(meta_elnr_s//'.CESL', 'L', jv_elnr_cesl)
    call jeveuo(meta_elnr_s//'.CESV', 'L', vr=v_elnr_cesv)
!
! - Fill META_ELNR_S on reduced list of internal variables
!
    do mapZoneNume = 1, mapNbZone
        nbCellInZone = comporInfoZone(mapZoneNume)
        if (nbCellInZone .ne. 0) then
! --------- Get object to link zone to internal variables
            call codent(mapZoneNume, 'G', saux08)
            vari_link = base_name//saux08
            call jeveuo(vari_link, 'L', vi=v_vari_link)

! --------- Access to current zone in CARTE
            affe_type = v_compor_desc(1+3+(mapZoneNume-1)*2)
            affe_indx = v_compor_desc(1+4+(mapZoneNume-1)*2)
            if (affe_type .eq. 3) then
                nb_elem = v_compor_lima_lc(1+affe_indx)-v_compor_lima_lc(affe_indx)
                posit = v_compor_lima_lc(affe_indx)
            elseif (affe_type .eq. 1) then
                nb_elem = nbCellMesh
                posit = 0
            else
                ASSERT(.false.)
            end if
            call jelira(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume), 'LONMAX', nb_vari_zone)
!
! --------- Loop on elements in zone of CARTE
!
            do i_elem = 1, nb_elem
                if (affe_type .eq. 3) then
                    nume_elem = v_compor_lima(posit+i_elem-1)
                elseif (affe_type .eq. 1) then
                    nume_elem = i_elem
                else
                    ASSERT(.false.)
                end if
                nb_pt = zi(jv_elno_cesd-1+5+4*(nume_elem-1)+1)
                nb_spt = zi(jv_elno_cesd-1+5+4*(nume_elem-1)+2)
                nb_vari = zi(jv_elno_cesd-1+5+4*(nume_elem-1)+3)
                do i_pt = 1, nb_pt
                    do i_spt = 1, nb_spt
                        do i_vari = 1, nb_vari
                            call cesexi('C', jv_elno_cesd, jv_elno_cesl, nume_elem, i_pt, &
                                        i_spt, i_vari, jv_elno)
                            if (jv_elno .gt. 0 .and. i_vari .le. nb_vari_zone) then
                                i_vari_redu = v_vari_link(i_vari)
                                if (i_vari_redu .ne. 0) then
                                    call cesexi('C', jv_elnr_cesd, jv_elnr_cesl, nume_elem, i_pt, &
                                                i_spt, i_vari_redu, jv_elnr)
                                    ASSERT(jv_elnr .ne. 0)
                                    jv_elnr = abs(jv_elnr)
                                    v_elnr_cesv(jv_elnr) = v_elga_cesv(jv_elno)
                                    zl(jv_elnr_cesl-1+jv_elnr) = .true.
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end if
    end do
!
! - Transform META_ELNR_S in META_ELNR
!
    nomres = field_med(1:8)//'META_ELNO_NOMME'
    call cescel(meta_elnr_s, ligrel, ' ', ' ', 'OUI', &
                nume_elem, 'V', meta_elnr, 'F', codret_dummy)
!
! - Write in MED file
!
    call irceme(ifi, nomres, meta_elnr, field_loca, ligrel, &
                nb_cmp_sele, cmp_name_sele, label_med, partie, numpt, &
                instan, nume_store, nbmaec, limaec, cara_elem, &
                cara_elem, field_type, nbCmpDyna, lfichUniq, codret)
!
99  continue

! - Cleaning
    call detrsd('CHAM_ELEM_S', meta_elno_s)
    call detrsd('CHAM_ELEM_S', meta_elnr_s)
    call detrsd('CHAM_ELEM_S', meta_elnr)
    call jedetr(comporInfo(1:19)//'.ZONE')
    call jedetr(comporInfo(1:19)//'.INFO')
    call jedetr(comporInfo(1:19)//'.ELEM')
    call jedetr(comporInfo(1:19)//'.RELA')
    call jedetc('V', comporInfo(1:19)//'.VARI', 1)
    call jedetc('V', comporInfo(1:19)//'.PHAS', 1)
    call jedetr(vari_redu)
    call jedetr(label_vxx)
    call jedetr(label_med)
    do mapZoneNume = 1, mapNbZone
        call codent(mapZoneNume, 'G', saux08)
        vari_link = base_name//saux08
        call jedetr(vari_link)
    end do
!
    call jedema()
!
end subroutine
