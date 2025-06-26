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
subroutine romModeSave(resultName, numeMode, &
                       fieldIden, mode, &
                       modeValeR_, modeValeC_, &
                       modeSing_, numeSlice_, &
                       nbSnap_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/jeveuo.h"
#include "asterfort/romModeParaSave.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
!
    character(len=8), intent(in) :: resultName
    integer(kind=8), intent(in) :: numeMode
    type(ROM_DS_Field), intent(in) :: mode
    character(len=24), intent(in) :: fieldIden
    real(kind=8), optional, pointer :: modeValeR_(:)
    complex(kind=8), optional, pointer :: modeValeC_(:)
    real(kind=8), optional, intent(in) :: modeSing_
    integer(kind=8), optional, intent(in) :: numeSlice_, nbSnap_
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Save mode
!
! --------------------------------------------------------------------------------------------------
!
! In  resultName       : name of result datastructure to save mode
! In  numeMode         : index of mode
! In  fieldIden        : identificator of field (name in results datastructure)
! In  mode             : mode to save
! In  modeValeR        : singular vector for mode (real)
! In  modeValeC        : singular vector for mode (complex)
! In  modeSing         : singular value for mode
! In  numeSlice        : index of slice (for lineic mode)
! In  nbSnap           : number of snapshots used to construct base
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, numeSlice, nbSnap, nbEqua
    real(kind=8) :: modeSing
    real(kind=8), pointer :: resuFieldValeR(:) => null()
    complex(kind=8), pointer :: resuFieldValeC(:) => null()
    character(len=8) :: model
    character(len=24) :: modeSymbName, resultField, fieldRefe
    character(len=4) :: fieldSupp
!
! --------------------------------------------------------------------------------------------------
!
    numeSlice = 0
    if (present(numeSlice_)) then
        numeSlice = numeSlice_
    end if
    modeSing = 0.d0
    if (present(modeSing_)) then
        modeSing = modeSing_
    end if
    nbSnap = 0
    if (present(nbSnap_)) then
        nbSnap = nbSnap_
    end if
!
! - Get parameters from mode
!
    nbEqua = mode%nbEqua
    modeSymbName = mode%fieldName
    fieldRefe = mode%fieldRefe
    fieldSupp = mode%fieldSupp
    model = mode%model
!
! - Get field from result to save mode
!
    call rsexch(' ', resultName, fieldIden, numeMode, resultField, iret)
    ASSERT(iret .eq. 100 .or. iret .eq. 0 .or. iret .eq. 110)
!
! - Extend results datastructure
!
    if (iret .eq. 110) then
        call rsagsd(resultName, 0)
        call rsexch(' ', resultName, fieldIden, numeMode, resultField, iret)
        ASSERT(iret .eq. 100 .or. iret .eq. 0)
    end if
!
! - Create a new field for mode
!
    if (iret .eq. 100) then
        call copisd('CHAMP_GD', 'G', fieldRefe, resultField)
    end if
!
! - Access to values of mode
!
    if (fieldSupp .eq. 'NOEU') then
        if (present(modeValeC_)) then
            call jeveuo(resultField(1:19)//'.VALE', 'E', vc=resuFieldValeC)
        else
            call jeveuo(resultField(1:19)//'.VALE', 'E', vr=resuFieldValeR)
        end if
    elseif (fieldSupp .eq. 'ELGA') then
        if (present(modeValeC_)) then
            call jeveuo(resultField(1:19)//'.CELV', 'E', vc=resuFieldValeC)
        else
            call jeveuo(resultField(1:19)//'.CELV', 'E', vr=resuFieldValeR)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Copy values of mode
!
    if (present(modeValeC_)) then
        resuFieldValeC(:) = modeValeC_(1:nbEqua)
    else
        resuFieldValeR(:) = modeValeR_(1:nbEqua)
    end if
    call rsnoch(resultName, fieldIden, numeMode)
!
! - Save parameters
!
    call romModeParaSave(resultName, numeMode, &
                         model, modeSymbName, &
                         modeSing, numeSlice, &
                         nbSnap)
!
end subroutine
