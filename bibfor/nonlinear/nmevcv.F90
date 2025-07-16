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
subroutine nmevcv(sderro, list_func_acti, loopName)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmeceb.h"
#include "asterfort/nmerge.h"
#include "asterfort/nmlecv.h"
#include "asterfort/NonLinear_type.h"
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=24), intent(in) :: sderro
    character(len=4), intent(in) :: loopName
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! EVALUATION DE LA CONVERGENCE
!
! --------------------------------------------------------------------------------------------------
!
! In  listFuncActi     : list of active functionnalities
! In  sderro           : name of datastructure for events in algorithm
! In  loopName         : name of loop
!               'RESI' - BOUCLE SUR LES RESIDUS D'EQUILIBRE
!               'NEWT' - BOUCLE DE NEWTON
!               'FIXE' - BOUCLE DE POINT FIXE
!               'INST' - BOUCLE SUR LES PAS DE TEMPS
!               'CALC' - CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: cveven
    integer(kind=8) :: iEvent
    character(len=24) :: eventENIVJv, eventENOMJv, eventEFCTJv
    character(len=9) :: eventName, eventLevel
    character(len=24) :: eventActiFunc
    aster_logical :: dv, cv, lfonc, withAnd, withOr
    character(len=4) :: loopState
    aster_logical :: convEqui
    aster_logical :: cvresi, cvnewt, cvfixe, cvinst
    character(len=16), pointer :: eventENIV(:) => null()
    character(len=16), pointer :: eventENOM(:) => null()
    character(len=24), pointer :: eventEFCT(:) => null()
    character(len=24) :: eventCONVJv
    integer(kind=8), pointer :: eventCONV(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to datastructure
    eventENOMJv = sderro(1:19)//'.ENOM'
    eventENIVJv = sderro(1:19)//'.ENIV'
    eventEFCTJv = sderro(1:19)//'.EFCT'
    eventCONVJv = sderro(1:19)//'.CONV'
    call jeveuo(eventENOMJv, 'L', vk16=eventENOM)
    call jeveuo(eventENIVJv, 'L', vk16=eventENIV)
    call jeveuo(eventEFCTJv, 'L', vk24=eventEFCT)
    call jeveuo(eventCONVJv, 'L', vi=eventCONV)

! - Evaluate convergence
    cveven = ASTER_TRUE
    convEqui = ASTER_TRUE
!
    if (eventCONV(NB_LOOP+1) .eq. 1) then
        withAnd = ASTER_FALSE
        withOr = ASTER_TRUE
    else
        withAnd = ASTER_TRUE
        withOr = ASTER_FALSE
    end if
!
    do iEvent = 1, ZEVEN
        eventLevel = eventENIV(iEvent) (1:9)
        eventName = eventENOM(iEvent) (1:9)
        eventActiFunc = eventEFCT(iEvent)
        dv = ASTER_FALSE
        cv = ASTER_TRUE
        if (eventLevel .eq. 'CONV_'//loopName) then
            if (eventActiFunc .eq. ' ') then
                lfonc = ASTER_TRUE
            else
                lfonc = isfonc(list_func_acti, eventActiFunc)
            end if
            if (eventName(1:4) .eq. 'DIVE') then
                call nmerge(sderro, eventName, dv)
                cv = ASTER_TRUE
            else if (eventName(1:4) .eq. 'CONV') then
                call nmerge(sderro, eventName, cv)
                cv = cv .and. lfonc
                dv = ASTER_FALSE
            end if

! --------- For law between residuals
            if (eventName .eq. 'DIVE_RELA' .or. &
                eventName .eq. 'DIVE_MAXI' .or. &
                eventName .eq. 'DIVE_REFE' .or. &
                eventName .eq. 'DIVE_COMP') then

                if (eventName .eq. 'DIVE_RELA') then
                    if (withAnd) then
                        convEqui = (.not. dv) .and. cv
                    else if (withOr) then
                        convEqui = (.not. dv) .and. lfonc .and. cv
                    end if
                else
                    if (withAnd) then
                        convEqui = convEqui .and. (.not. dv) .and. cv
                    else if (withOr) then
                        convEqui = convEqui .or. ((.not. dv) .and. lfonc .and. cv)
                    end if
                end if
            else
                cveven = cveven .and. (.not. dv) .and. cv
            end if
        end if

    end do
    cveven = cveven .and. convEqui

! - Get state of previous loops
    call nmlecv(sderro, 'RESI', cvresi)
    call nmlecv(sderro, 'NEWT', cvnewt)
    call nmlecv(sderro, 'FIXE', cvfixe)
    call nmlecv(sderro, 'INST', cvinst)
    if (loopName .eq. 'NEWT') then
        cveven = cveven .and. cvresi
    end if
    if (loopName .eq. 'FIXE') then
        cveven = cveven .and. cvnewt .and. cvresi
    end if
    if (loopName .eq. 'INST') then
        cveven = cveven .and. cvfixe .and. cvnewt .and. cvresi
    end if
    if (loopName .eq. 'CALC') then
        cveven = cveven .and. cvinst .and. cvfixe .and. cvnewt .and. cvresi
    end if

! - Set state of loop (continue loop/stop loop because converges)
    loopState = 'CONT'
    if (cveven) then
        loopState = 'CONV'
    end if

! - Save state of loop
    call nmeceb(sderro, loopName, loopState)
!
    call jedema()
end subroutine
