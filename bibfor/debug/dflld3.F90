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
subroutine dflld3(sdlist, iAdap)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: sdlist
    integer(kind=8), intent(in) :: iAdap
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_LIST_INST
!
! Print for next time step evaluation
!
! --------------------------------------------------------------------------------------------------
!
! In  sdlist           : name of DEFI_LIST_INST datastructure
! In  iAdap            : current scheme of adaptation
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdlistAEvenrName
    real(kind=8), pointer :: sdlistAEvenr(:) => null()
    character(len=24) :: sdlistATplurName
    real(kind=8), pointer :: sdlistATplur(:) => null()
    character(len=24) :: sdlistATplukName
    character(len=16), pointer :: sdlistATpluk(:) => null()
    integer(kind=8) :: action_type, nb_iter_newton_ref
    real(kind=8) :: pcent_augm, vale_ref
    character(len=16):: nom_cham, nom_cmp, crit_cmp
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to datastructures
    sdlistAEvenrName = sdlist(1:8)//'.ADAP.EVENR'
    sdlistATplurName = sdlist(1:8)//'.ADAP.TPLUR'
    sdlistATplukName = sdlist(1:8)//'.ADAP.TPLUK'
    call jeveuo(sdlistAEvenrName, 'L', vr=sdlistAEvenr)
    call jeveuo(sdlistATplurName, 'L', vr=sdlistATplur)
    call jeveuo(sdlistATplukName, 'L', vk16=sdlistATpluk)

! - Print
    action_type = nint(sdlistATplur(SIZE_LATPR*(iAdap-1)+1))
    pcent_augm = sdlistATplur(SIZE_LATPR*(iAdap-1)+2)
    vale_ref = sdlistATplur(SIZE_LATPR*(iAdap-1)+3)
    nom_cham = sdlistATpluk(SIZE_LATPK*(iAdap-1)+2)
    nom_cmp = sdlistATpluk(SIZE_LATPK*(iAdap-1)+3)
    crit_cmp = 'GE'
    nb_iter_newton_ref = nint(sdlistATplur(SIZE_LATPR*(iAdap-1)+5))
    if (action_type .eq. ADAP_ACT_FIXE) then
        call utmess('I', 'DISCRETISATION3_80')
        call utmess('I', 'DISCRETISATION3_84', sr=pcent_augm)
    elseif (action_type .eq. ADAP_ACT_INCR_QUANT) then
        call utmess('I', 'DISCRETISATION3_81')
        call utmess('I', 'DISCRETISATION3_21', &
                    nk=3, valk=[nom_cham, nom_cmp, crit_cmp], &
                    sr=vale_ref)
    elseif (action_type .eq. ADAP_ACT_ITER) then
        call utmess('I', 'DISCRETISATION3_82')
        call utmess('I', 'DISCRETISATION3_85', si=nb_iter_newton_ref)
    elseif (action_type .eq. ADAP_ACT_IMPLEX) then
        call utmess('I', 'DISCRETISATION3_83')
    end if
!
    call jedema()
end subroutine
