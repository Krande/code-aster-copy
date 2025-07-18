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

subroutine apstoc(ds_contact, nb_pair, list_pair, list_nbptit, list_ptitsl, &
                  list_ptitma, list_ptgama)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/wkvect.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/as_deallocate.h"
!
!
    type(NL_DS_Contact), intent(inout) :: ds_contact
    integer(kind=8), intent(in):: nb_pair
    integer(kind=8), pointer :: list_pair(:)
    integer(kind=8), pointer :: list_nbptit(:)
    real(kind=8), pointer :: list_ptitsl(:)
    real(kind=8), pointer :: list_ptitma(:)
    real(kind=8), pointer :: list_ptgama(:)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Save pairing information in sdappa data structure
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_contact       : datastructure for contact management
! In  nb_pair          : number of pairs in contact zone
! IO  list_pair        : list of pairs in contact zone
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: sdappa
    integer(kind=8), pointer :: v_sdappa_apli(:) => null()
    integer(kind=8), pointer :: v_sdappa_apli2(:) => null()
    real(kind=8), pointer :: v_sdappa_apli3(:) => null()
    real(kind=8), pointer :: v_sdappa_apli4(:) => null()
    real(kind=8), pointer :: v_sdappa_apli5(:) => null()
    character(len=24) :: sdappa_apli, sdappa_apli2, sdappa_apli3, sdappa_apli4, sdappa_apli5
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    ds_contact%nb_cont_pair = nb_pair
    if (nb_pair .ne. 0) then
        sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
        sdappa_apli = sdappa(1:19)//'.APLI'
        sdappa_apli2 = sdappa(1:19)//'.APNP'
        sdappa_apli3 = sdappa(1:19)//'.APTS'
        sdappa_apli4 = sdappa(1:19)//'.APTM'
        sdappa_apli5 = sdappa(1:19)//'.AP2M'
        call jedetr(sdappa_apli)
        call jedetr(sdappa_apli2)
        call jedetr(sdappa_apli3)
        call jedetr(sdappa_apli4)
        call jedetr(sdappa_apli5)
        call wkvect(sdappa_apli, 'V V I', 3*nb_pair, vi=v_sdappa_apli)
        call wkvect(sdappa_apli2, 'V V I', nb_pair, vi=v_sdappa_apli2)
        call wkvect(sdappa_apli3, 'V V R', 16*nb_pair, vr=v_sdappa_apli3)
        call wkvect(sdappa_apli4, 'V V R', 16*nb_pair, vr=v_sdappa_apli4)
        call wkvect(sdappa_apli5, 'V V R', 72*nb_pair, vr=v_sdappa_apli5)
        v_sdappa_apli(1:3*nb_pair) = list_pair(1:3*nb_pair)
        v_sdappa_apli2(1:nb_pair) = list_nbptit(1:nb_pair)
        v_sdappa_apli3(1:16*nb_pair) = list_ptitsl(1:16*nb_pair)
        v_sdappa_apli4(1:16*nb_pair) = list_ptitma(1:16*nb_pair)
        v_sdappa_apli5(1:72*nb_pair) = list_ptgama(1:72*nb_pair)
    else
        sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
        sdappa_apli = sdappa(1:19)//'.APLI'
        sdappa_apli2 = sdappa(1:19)//'.APNP'
        sdappa_apli3 = sdappa(1:19)//'.APTS'
        sdappa_apli4 = sdappa(1:19)//'.APTM'
        sdappa_apli5 = sdappa(1:19)//'.AP2M'
        call jedetr(sdappa_apli)
        call jedetr(sdappa_apli2)
        call jedetr(sdappa_apli3)
        call jedetr(sdappa_apli4)
        call jedetr(sdappa_apli5)
    end if
    AS_DEALLOCATE(vi=list_pair)
    AS_DEALLOCATE(vi=list_nbptit)
    AS_DEALLOCATE(vr=list_ptitsl)
    AS_DEALLOCATE(vr=list_ptitma)
    AS_DEALLOCATE(vr=list_ptgama)

    call jedema()
end subroutine
