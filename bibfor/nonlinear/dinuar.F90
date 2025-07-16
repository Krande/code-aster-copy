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

subroutine dinuar(result, sddisc, numeInst, force, &
                  numeStoring, nume_reuse_, lStoringInitState_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/diinst.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcrpo.h"
#include "asterfort/rsadpa.h"
!
    character(len=8), intent(in) :: result
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeInst
    aster_logical, intent(in) :: force
    integer(kind=8), intent(out) :: numeStoring
    integer(kind=8), optional, intent(out) :: nume_reuse_
    aster_logical, intent(in), optional :: lStoringInitState_
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Get storing index
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of datastructure for results
! In  sddisc           : datastructure for time discretization
! In  numeInst         : index of current time step
! In  force            : to "force" storing (ex.: error)
! Out numeStoring      : index to store in results
! Out numeReuseCalc    : index for reuse rsults datastructure
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdarchAinfJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    integer(kind=8) :: numeReuseCalc, jv_para
    real(kind=8) :: time_curr, time_prev
    aster_logical :: l_store, lStoringInitState
    character(len=19) :: sdarch
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    l_store = ASTER_FALSE
    lStoringInitState = ASTER_FALSE
    if (present(lStoringInitState_)) then
        lStoringInitState = lStoringInitState_
    end if
    numeStoring = -1
    numeReuseCalc = -1

! - Acces to storing objects
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    call jeveuo(sdarchAinfJv, 'E', vi=sdarchAinf)

! - Initial storing => already save
    if (numeInst .eq. 0) then
        l_store = ASTER_TRUE
    end if

! - Other time step => to save or not ?
    if (numeInst .ne. 0) then
        time_curr = diinst(sddisc, numeInst)
        call nmcrpo(sdarch, numeInst, time_curr, l_store)
    end if

! - "forced" storing
    if (force) then
        l_store = .true.
    end if

! - Storing index
    if (l_store) then
        numeStoring = sdarchAinf(1)
    else
        numeStoring = -1
    end if

! - REUSE for PARA_CALC table
    numeReuseCalc = sdarchAinf(3)

! - Already stored ?
    if (numeStoring .ge. 2) then
        call rsadpa(result, 'L', 1, 'INST', numeStoring-1, 0, sjv=jv_para)
        time_prev = zr(jv_para)
        if (time_curr .le. time_prev) then
            numeStoring = -1
            l_store = .false._1
        end if
    end if

! - Increase storing index
    if (l_store) then
        sdarchAinf(1) = sdarchAinf(1)+1
        if (lStoringInitState) then
            sdarchAinf(4) = 0
        else
            sdarchAinf(4) = sdarchAinf(4)+1
        end if
    end if
!
    if (present(nume_reuse_)) then
        nume_reuse_ = numeReuseCalc
    end if
!
    call jedema()
!
end subroutine
