! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine nmcrli(inst_init, list_inst, sddisc)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
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
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: list_inst
    real(kind=8), intent(in) :: inst_init
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Datastructures
!
! Time discretization datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  inst_init        : initial time if ETAT_INIT
! In  list_inst        : list of times from INCREMENT/LIST_INST
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    integer :: nume_ini, nume_end, nume_inst
    integer :: nb_inst_new, nb_inst, nbret
    real(kind=8) :: tole
    real(kind=8) :: dtmin, dt0
    aster_logical :: l_init_noexist
    character(len=24) :: list_inst_info
    character(len=24) :: list_inst_ditr
    character(len=16) :: list_inst_type, keywf
    character(len=24) :: sddisc_bcle
    integer, pointer :: v_sddisc_bcle(:) => null()
    character(len=24) :: sddisc_epil
    integer, pointer :: v_sddisc_epil(:) => null()
    character(len=19) :: list_inst_work
    real(kind=8), pointer :: v_list_work(:) => null()
    character(len=24) :: sddisc_dini
    integer, pointer :: v_sddisc_dini(:) => null()
    character(len=24) :: sddisc_iter
    integer, pointer :: v_sddisc_iter(:) => null()
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
!
! - Initializations
!
    l_init_noexist = .false.
    keywf = 'INCREMENT'
    list_inst_work = '&&NMCRLI.PROVLI'
!
! - Create loops object
! --- 1 - Newton (ITERAT)
! --- 2 - Time stepping (NUME_INST)
! --- 3 - Fixed loops (NIVEAU)
!
    sddisc_bcle = sddisc(1:19)//'.BCLE'
    call wkvect(sddisc_bcle, 'V V I', 3, vi=v_sddisc_bcle)
!
! - Create continuation choice object
!
    sddisc_epil = sddisc(1:19)//'.EPIL'
    call wkvect(sddisc_epil, 'V V I', 2, vi=v_sddisc_epil)
    v_sddisc_epil(1) = 1
    v_sddisc_epil(2) = 1
!
! - Type of list_inst
!
    call gettco(list_inst, list_inst_type)
    ASSERT(list_inst_type .ne. ' ')
!
! - Create list of times and information vector
!
    if (list_inst_type .eq. 'LISTR8_SDASTER') then
        call nmcrlm(list_inst, sddisc, list_inst_work)
    else if (list_inst_type .eq. 'LIST_INST') then
        sddisc_linf = sddisc(1:19)//'.LINF'
        list_inst_info = list_inst(1:8)//'.LIST.INFOR'
        list_inst_ditr = list_inst(1:8)//'.LIST.DITR'
        call jedup1(list_inst_ditr, 'V', list_inst_work)
        call jedup1(list_inst_info, 'V', sddisc_linf)
    end if
!
! - Get parameters
!
    call utdidt('L', sddisc, 'LIST', 'DTMIN', &
                valr_=dtmin)
    call utdidt('L', sddisc, 'LIST', 'NBINST', &
                vali_=nb_inst)
!
! - Acces to list of times
!
    call jeveuo(list_inst_work, 'L', vr=v_list_work)
!
! - Get parameters
!
    call getvr8(keywf, 'PRECISION', iocc=1, scal=tole, nbret=nbret)
    if (nbret == 0) then
        tole = 1d-6
    end if
    tole = abs(dtmin)*tole
!
! - Index of initial time
!
    call nmdini(keywf, list_inst_work, tole, &
                nb_inst, l_init_noexist, nume_ini)
!
! - Index of final time
!
    call nmdifi(keywf, list_inst_work, tole, nb_inst, nume_end)
!
! - Check
!
    if (nume_ini .ge. nume_end) then
        call utmess('F', 'DISCRETISATION_92')
    end if
!
! - Resize list of times
!
    call nmcrls(sddisc, list_inst_work, nume_ini, nume_end, l_init_noexist, &
                inst_init, nb_inst_new, dtmin)
!
! - Create object for subdividing time steps
!
    sddisc_dini = sddisc(1:19)//'.DINI'
    call wkvect(sddisc_dini, 'V V I', nb_inst_new, vi=v_sddisc_dini)
    do nume_inst = 1, nb_inst_new
        v_sddisc_dini(nume_inst) = 1
    end do
!
! - Create object for number of iterations
!
    sddisc_iter = sddisc(1:19)//'.ITER'
    call wkvect(sddisc_iter, 'V V I', nb_inst_new, vi=v_sddisc_iter)
!
! - Save parameters
!
    dt0 = diinst(sddisc, 1)-diinst(sddisc, 0)
    call utdidt('E', sddisc, 'LIST', 'DT-', &
                valr_=dt0)
    call utdidt('E', sddisc, 'LIST', 'NBINST', &
                vali_=nb_inst_new)
    call utdidt('E', sddisc, 'LIST', 'DTMIN', &
                valr_=dtmin)
!
! - Save object of time steps
!
    sddisc_ditr = sddisc(1:19)//'.DITR'
    sddisc_lipo = sddisc(1:19)//'.LIPO'
    call jedupo(sddisc_ditr, 'V', sddisc_lipo, .false._1)
!
! - Clean
!
    call jedetr(list_inst_work)
!
end subroutine
