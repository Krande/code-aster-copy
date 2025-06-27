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
subroutine nmarch(numeInst, model, ds_material, caraElem, listFuncActi, &
                  ds_print, sddisc, sdcrit, &
                  ds_measure, sderro, sddyna, sdpilo, ds_energy, &
                  ds_inout, ds_errorindic, ds_algorom_, lStoringInitState_)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use HHO_postpro_module, only: hhoPostDeplMeca
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/isfonc.h"
#include "asterfort/diinst.h"
#include "asterfort/dinuar.h"
#include "asterfort/nmarc0.h"
#include "asterfort/nmarce.h"
#include "asterfort/nmarpc.h"
#include "asterfort/nmcrpc.h"
#include "asterfort/nmfinp.h"
#include "asterfort/nmleeb.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpg.h"
#include "asterfort/romAlgoNLTableSave.h"
!
    integer(kind=8) :: listFuncActi(*)
    integer(kind=8) :: numeInst
    type(NL_DS_Print), intent(in) :: ds_print
    type(NL_DS_InOut), intent(in) :: ds_inout
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Energy), intent(in) :: ds_energy
    character(len=19) :: sddisc, sdcrit, sddyna, sdpilo
    character(len=24) :: sderro
    type(NL_DS_ErrorIndic), intent(in) :: ds_errorindic
    character(len=24) :: model, caraElem
    type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
    aster_logical, intent(in), optional :: lStoringInitState_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Storing results
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_inout         : datastructure for input/output management
! In  ds_print         : datastructure for printing parameters
! IN  NUMINS : NUMERO DE L'INSTANT
! IN  MODELE : NOM DU MODELEE
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  SDCRIT : VALEUR DES CRITERES DE CONVERGENCE
! In  ds_errorindic    : datastructure for error indicator
! IN  SDERRO : SD ERREUR
! IN  SDDYNA : SD DEDIEE A LA DYNAMIQUE
! IN  SDPILO : SD PILOTAGE
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_energy        : datastructure for energy management
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! In  ds_algorom       : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, numeStoring
    real(kind=8) :: timeCurr
    character(len=8) :: result
    aster_logical :: forceStoring, lprint, l_hho, lStoringInitState, lastTimeStep
    character(len=19) :: k19bid
    character(len=24) :: listLoadResu
    character(len=4) :: etcalc
    integer(kind=8) :: numeReuse
!
! --------------------------------------------------------------------------------------------------
!
    lStoringInitState = ASTER_FALSE
    if (present(lStoringInitState_)) then
        lStoringInitState = lStoringInitState_
    end if
    result = ds_inout%result
    listLoadResu = ds_inout%listLoadResu
    l_hho = isfonc(listFuncActi, 'HHO')

! - Loop state
    call nmleeb(sderro, 'CALC', etcalc)

! - Last step => storing
    forceStoring = ASTER_FALSE
    call nmfinp(sddisc, numeInst, lastTimeStep)
    forceStoring = lastTimeStep

! - Storing
    if (etcalc .eq. 'CONV') then
        forceStoring = .true.
    end if
    if (etcalc .eq. 'STOP') then
        forceStoring = .true.
    end if

! - Print timer
    call uttcpg('IMPR', 'INCR')

! - Get index for storing
    call dinuar(result, sddisc, numeInst, forceStoring, &
                numeStoring, numeReuse, lStoringInitState)

! - Current time
    timeCurr = diinst(sddisc, numeInst)

! - Save energy parameters in output table
    if (isfonc(listFuncActi, 'ENERGIE')) then
        call nmarpc(ds_energy, numeReuse, timeCurr)
    else
        call nmcrpc(ds_inout, numeReuse, timeCurr)
    end if

! - Print or not ?
    lprint = ds_print%l_print

! - Storing
    if (numeStoring .ge. 0) then
! ----- Begin timer
        call nmtime(ds_measure, 'Launch', 'Store')

! ----- Print head
        if (lprint) then
            call utmess('I', 'ARCHIVAGE_5')
        end if

! ----- Increased result datastructure if necessary
        call rsexch(' ', result, 'DEPL', numeStoring, k19bid, iret)
        if (iret .eq. 110) then
            call rsagsd(result, 0)
        end if

! ----- Storing parameters
        call nmarc0(result, model, ds_material, caraElem, listFuncActi, &
                    sdcrit, sddyna, ds_errorindic, &
                    sdpilo, listLoadResu, numeStoring, timeCurr)

! ----- Storing fields
        call nmarce(ds_inout, result, sddisc, timeCurr, numeStoring, &
                    forceStoring, ds_print)

! ----- If HHO, we compute a post_processing
        if (l_hho) then
            call hhoPostDeplMeca(model, result, numeStoring)
        end if

! ----- Storing reduced parameters table (ROM)
        if (present(ds_algorom_)) then
            if (ds_algorom_%l_rom) then
                if (numeStoring .gt. 0) then
                    call romAlgoNLTableSave(numeStoring, timeCurr, ds_algorom_)
                end if
            end if
        end if

! ----- End timer
        call nmtime(ds_measure, 'Stop', 'Store')
        call nmrinc(ds_measure, 'Store')
    end if
!
end subroutine
