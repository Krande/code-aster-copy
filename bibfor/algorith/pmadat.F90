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
subroutine pmadat(sddisc, numeInst, nbIter)
!
    implicit none
!
#include "asterc/r8maem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/compr8.h"
#include "asterfort/diadap.h"
#include "asterfort/diinst.h"
#include "asterfort/getAdapAction.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/nmdcei.h"
#include "asterfort/nmjalo.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "event_def.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeInst, nbIter
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! GESTION DE L'ADAPTATION DE PAS DE TEMPS
!   CALCUL DE DTPLUS
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! IN  NUMINS : NUMERO D'INSTANT
! IN  NBITER : NOMBRE D'ITERATIONS DE NEWTON
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: dtPlusJv = '&&NMADAP.DTPLUS'
    real(kind=8), pointer :: dtPlus(:) => null()
    real(kind=8), parameter :: prec = 1.d-12
    integer(kind=8), parameter :: inspas = 1
    integer(kind=8) :: nbAdap, iAdap
    character(len=19) :: listInstMethod
    real(kind=8) :: r8bid, dt, min, pasMini, pasMaxi, dtm, jalon
    real(kind=8) :: newdt, deltac
    real(kind=8) :: timeCurr, timeNew
    real(kind=8) :: insfin, insref
    integer(kind=8) :: nbInst, nbPasMaxi, actionType
    aster_logical :: ladap, uncrok, lOutOFTime
    character(len=24) :: tpsIterJv
    integer(kind=8), pointer :: tpsIter(:) => null()
    character(len=24) :: tpsAextJv
    real(kind=8), pointer :: tpsAext(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Get parameters
    call utdidt('L', sddisc, 'LIST', 'PAS_MINI', valr_=pasMini)
    call utdidt('L', sddisc, 'LIST', 'PAS_MAXI', valr_=pasMaxi)
    call utdidt('L', sddisc, 'LIST', 'METHODE', valk_=listInstMethod)
    call utdidt('L', sddisc, 'LIST', 'NBINST', vali_=nbInst)
    call utdidt('L', sddisc, 'LIST', 'NB_PAS_MAXI', vali_=nbPasMaxi)
    tpsIterJv = sddisc(1:19)//'.ITER'

    if (listInstMethod .ne. 'MANUEL') then
! ----- Current time
        timeCurr = diinst(sddisc, numeInst)

! ----- Seek for next time step
        call nmjalo(sddisc, timeCurr, prec, jalon, lOutOFTime)
        if (.not. lOutOFTime) then
            call utdidt('L', sddisc, 'LIST', 'NADAPT', vali_=nbAdap)
            call wkvect(dtPlusJv, 'V V R', nbAdap, vr=dtPlus)
            call utdidt('L', sddisc, 'LIST', 'DT-', valr_=dtm)
            call jeveuo(tpsIterJv, 'L', vi=tpsIter)
            tpsIter(numeInst) = nbIter
            call juveca(tpsIterJv, nbInst+1)

            call utmess('I', 'ADAPTATION_1')
            do iAdap = 1, nbAdap
                call getAdapAction(sddisc, iAdap, actionType)
                dtPlus(iAdap) = r8vide()
                ladap = diadap(sddisc, iAdap)
                if (ladap) then
                    call utmess("F", "COMPOR2_27")
                end if
                newdt = dtPlus(iAdap)
                if (newdt .ne. r8vide()) then
                    call utmess('I', 'ADAPTATION_2', sk=adapActionKeyword(actionType), sr=newdt)
                else
                    call utmess('I', 'ADAPTATION_3', sk=adapActionKeyword(actionType))
                end if
            end do

! --------- Get minimal DeltaT
            dt = r8maem()
            uncrok = ASTER_FALSE
            do iAdap = 1, nbAdap
                newdt = dtPlus(iAdap)
                if (newdt .ne. r8vide()) then
                    dt = min(dt, newdt)
                    uncrok = ASTER_TRUE
                end if
            end do
            if (uncrok) then
                call utmess('I', 'ADAPTATION_5', sr=dt)
            else
                dt = dtm
                call utmess('I', 'ADAPTATION_4', sr=dt)
            end if

            if (dt .gt. pasMaxi) then
                call utmess('I', 'ADAPTATION_12', nr=2, valr=[dt, pasMaxi])
                dt = pasMaxi
            end if

            tpsAextJv = sddisc(1:19)//'.AEXT'
            call jeveuo(tpsAextJv, 'E', vr=tpsAext)
            insref = tpsAext(1)
            deltac = tpsAext(2)
            insfin = tpsAext(3)
            if (insref .ne. r8vide()) then
                if (timeCurr .le. insfin) then
                    dt = deltac
                    call utmess('I', 'ADAPTATION_10', sr=dt)
                else
                    tpsAext(1) = r8vide()
                end if
            end if

            if (compr8(timeCurr+dt, 'GT', jalon, prec, 1)) then
                dt = jalon-timeCurr
            else if (compr8(timeCurr+dt, 'GT', jalon-pasMini, prec, 1)) then
                dt = jalon-timeCurr
                call utdidt('E', sddisc, 'LIST', 'DT-', valr_=dt)
            else
                call utdidt('E', sddisc, 'LIST', 'DT-', valr_=dt)
            end if
            call utmess('I', 'ADAPTATION_6', sr=dt)
            if (dt .lt. pasMini) then
                call utmess('F', 'ADAPTATION_11', sr=dt)
            end if

            timeNew = timeCurr+dt
            call nmdcei(sddisc, numeInst, [timeNew], nbInst, inspas, 'ADAP', r8bid)
        end if
    end if
!
    call jedetr(dtPlusJv)
    call jedema()
end subroutine
