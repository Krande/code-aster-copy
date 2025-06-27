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
subroutine nmcrli(listInst, sddisc)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/gettco.h"
#include "asterfort/getvr8.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedup1.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcrlm.h"
#include "asterfort/nmcrls.h"
#include "asterfort/nmdifi.h"
#include "asterfort/nmdini.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=19), intent(in) :: sddisc, listInst
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Datastructures
!
! Time discretization datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc          : datastructure for time discretization
! In  listInst        : list of times from INCREMENT/LIST_INST
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'INCREMENT'
    character(len=19), parameter :: listInstWorkJv = '&&NMCRLI.PROVLI'
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: numeInit, numeEnd
    integer(kind=8) :: nbInstNew, nbInst, nbret
    real(kind=8) :: tole, dtmin, dt0, instInit
    character(len=24) :: list_inst_info, list_inst_ditr
    character(len=16) :: list_inst_type
    character(len=24) :: sddisc_bcle
    integer(kind=8), pointer :: v_sddisc_bcle(:) => null()
    real(kind=8), pointer :: listInstWork(:) => null()
    character(len=24) :: sddisc_dini
    integer(kind=8), pointer :: v_sddisc_dini(:) => null()
    character(len=24) :: sddisc_iter
    integer(kind=8), pointer :: v_sddisc_iter(:) => null()
    character(len=24) :: sddisc_lipo
    character(len=24) :: sddisc_ditr
    character(len=24) :: sddisc_linf
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_15')
    end if

! - Create loops object
! --- 1 - Newton (ITERAT)
! --- 2 - Time stepping (NUME_INST)
! --- 3 - Fixed loops (NIVEAU)
    sddisc_bcle = sddisc(1:19)//'.BCLE'
    call wkvect(sddisc_bcle, 'V V I', 3, vi=v_sddisc_bcle)

! - Type of listInst
    call gettco(listInst, list_inst_type)
    ASSERT(list_inst_type .ne. ' ')

! - Create list of times and information vector
    if (list_inst_type .eq. 'LISTR8_SDASTER') then
        call nmcrlm(listInst, sddisc, listInstWorkJv)
    else if (list_inst_type .eq. 'LIST_INST') then
        sddisc_linf = sddisc(1:19)//'.LINF'
        list_inst_info = listInst(1:8)//'.LIST.INFOR'
        list_inst_ditr = listInst(1:8)//'.LIST.DITR'
        call jedup1(list_inst_ditr, 'V', listInstWorkJv)
        call jedup1(list_inst_info, 'V', sddisc_linf)
    end if

! - Get parameters
    call utdidt('L', sddisc, 'LIST', 'DTMIN', valr_=dtmin)
    call utdidt('L', sddisc, 'LIST', 'NBINST', vali_=nbInst)

! - Acces to list of times
    call jeveuo(listInstWorkJv, 'L', vr=listInstWork)

! - Get parameters
    call getvr8(factorKeyword, 'PRECISION', iocc=1, scal=tole, nbret=nbret)
    if (nbret == 0) then
        tole = 1d-6
    end if
    tole = abs(dtmin)*tole

! - Index of initial time
    call nmdini(factorKeyword, listInstWorkJv, tole, &
                nbInst, numeInit, instInit)

! - Index of final time
    call nmdifi(factorKeyword, listInstWorkJv, tole, nbInst, numeEnd)

! - Check
    if (numeInit .ge. numeEnd) then
        call utmess('F', 'DISCRETISATION_92')
    end if

! - Resize list of times
    call nmcrls(sddisc, listInstWorkJv, numeInit, numeEnd, &
                nbInstNew, dtmin)

! - Create object for subdividing time steps
    sddisc_dini = sddisc(1:19)//'.DINI'
    call wkvect(sddisc_dini, 'V V I', nbInstNew, vi=v_sddisc_dini)
    v_sddisc_dini(1:nbInstNew) = 1

! - Create object for number of iterations
    sddisc_iter = sddisc(1:19)//'.ITER'
    call wkvect(sddisc_iter, 'V V I', nbInstNew, vi=v_sddisc_iter)

! - Save parameters
    dt0 = diinst(sddisc, 1)-diinst(sddisc, 0)
    call utdidt('E', sddisc, 'LIST', 'DT-', valr_=dt0)
    call utdidt('E', sddisc, 'LIST', 'NBINST', vali_=nbInstNew)
    call utdidt('E', sddisc, 'LIST', 'DTMIN', valr_=dtmin)

! - Save object of time steps
    sddisc_ditr = sddisc(1:19)//'.DITR'
    sddisc_lipo = sddisc(1:19)//'.LIPO'
    call jedupo(sddisc_ditr, 'V', sddisc_lipo, .false._1)

! - Clean
    call jedetr(listInstWorkJv)
!
end subroutine
