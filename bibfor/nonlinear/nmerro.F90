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
subroutine nmerro(sderro, ds_measure, nume_inst)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterc/etausr.h"
#include "asterfort/nmecev.h"
#include "asterfort/nmerge.h"
#include "asterfort/sigusr.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: sderro
    type(NL_DS_Measure), intent(in) :: ds_measure
    integer(kind=8), intent(in) :: nume_inst
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Error management
!
! Write messages for errors
!
! --------------------------------------------------------------------------------------------------
!
! In  sderro           : datastructure for errors during algorithm
! In  ds_measure       : datastructure for measure and statistics management
! In  nume_inst        : index of current time step
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: rtab(2)
    integer(kind=8) :: itab(2)
    aster_logical :: echldc, echeq1, echeq2, echco1, echco2, echpil
    aster_logical :: mtcpui, mtcpup, itemax
    aster_logical :: echpfg, echpff, echpfc
    aster_logical :: errres, err_appa
    real(kind=8) :: remain_time, iter_mean_time, step_mean_time
    character(len=16) :: valk(2)
    integer(kind=8) :: failType, actionType
!
! --------------------------------------------------------------------------------------------------
!
    if (etausr() .eq. 1) then
        call sigusr()
    end if

! - Get times
    remain_time = ds_measure%step_remain_time
    iter_mean_time = ds_measure%iter_mean_time
    step_mean_time = ds_measure%step_mean_time

! --- RECUPERE LES CODES ERREURS ACTIFS
    call nmerge(sderro, 'ERRE_INTE', echldc)
    call nmerge(sderro, 'ERRE_PILO', echpil)
    call nmerge(sderro, 'ERRE_FACS', echeq1)
    call nmerge(sderro, 'ERRE_FACT', echeq2)
    call nmerge(sderro, 'ERRE_CTD1', echco1)
    call nmerge(sderro, 'ERRE_CTD2', echco2)
    call nmerge(sderro, 'ERRE_TIMN', mtcpui)
    call nmerge(sderro, 'ERRE_TIMP', mtcpup)
    call nmerge(sderro, 'ITER_MAXI', itemax)
    call nmerge(sderro, 'ERRE_CTCG', echpfg)
    call nmerge(sderro, 'ERRE_CTCF', echpff)
    call nmerge(sderro, 'ERRE_CTCC', echpfc)
    call nmerge(sderro, 'SOLV_ITMX', errres)
    call nmerge(sderro, 'ERRE_APPA', err_appa)

! --- LANCEE EXCEPTIONS
    if (mtcpui) then
        itab(1) = nume_inst
        rtab(1) = iter_mean_time
        rtab(2) = remain_time
        call utmess('Z', 'MECANONLINE9_1', si=itab(1), nr=2, valr=rtab, &
                    num_except=ASTER_TIMELIMIT_ERROR)
    else if (mtcpup) then
        itab(1) = nume_inst
        rtab(1) = step_mean_time
        rtab(2) = remain_time
        call utmess('Z', 'MECANONLINE9_2', si=itab(1), nr=2, valr=rtab, &
                    num_except=ASTER_TIMELIMIT_ERROR)
    else if (echldc) then
        call utmess('Z', 'MECANONLINE9_3', num_except=ASTER_INTEGRATION_ERROR)
    else if (echeq1 .or. echeq2) then
        call utmess('Z', 'MECANONLINE9_4', num_except=ASTER_SOLVER_ERROR)
    else if (echco1) then
        call utmess('Z', 'MECANONLINE9_5', num_except=ASTER_CONTACT_ERROR)
    else if (echco2) then
        call utmess('Z', 'MECANONLINE9_6', num_except=ASTER_SOLVER_ERROR)
    else if (itemax) then
        call utmess('Z', 'MECANONLINE9_7', num_except=ASTER_CONVERGENCE_ERROR)
    else if (echpil) then
        call utmess('Z', 'MECANONLINE9_8', num_except=ASTER_CONVERGENCE_ERROR)
    else if (echpfg) then
        call utmess('Z', 'MECANONLINE9_9', num_except=ASTER_CONTACT_ERROR)
    else if (echpff) then
        call utmess('Z', 'MECANONLINE9_10', num_except=ASTER_CONTACT_ERROR)
    else if (echpfc) then
        call utmess('Z', 'MECANONLINE9_11', num_except=ASTER_CONTACT_ERROR)
    else if (errres) then
        call utmess('Z', 'MECANONLINE9_12', num_except=ASTER_SOLVER_ERROR)
    else if (err_appa) then
        call utmess('Z', 'MECANONLINE9_13', num_except=ASTER_CONTACT_ERROR)
    else
        call nmecev(sderro, 'L', failType, actionType)
        valk(1) = failActionKeyword(actionType)
        valk(2) = failEventKeyword(failType)
        if (actionType .eq. FAIL_ACT_STOP) then
            call utmess('Z', 'MECANONLINE9_51', sk=failEventKeyword(failType), &
                        num_except=ASTER_CONVERGENCE_ERROR)
        else
            call utmess('Z', 'MECANONLINE9_50', nk=2, valk=valk, &
                        num_except=ASTER_CONVERGENCE_ERROR)
        end if
    end if
!
end subroutine
