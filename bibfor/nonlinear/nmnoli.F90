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
subroutine nmnoli(sddisc, sderro, ds_print, sdcrit, &
                  fonact, sddyna, modele, ds_material, &
                  carele, sdpilo, ds_measure, ds_energy, ds_inout, &
                  ds_errorindic)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/isfonc.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmarch.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsrusd.h"
#include "asterfort/utmess.h"
!
    character(len=19) :: sddisc, sdcrit, sddyna, sdpilo
    type(NL_DS_Energy), intent(in) :: ds_energy
    character(len=24) :: sderro
    character(len=24) :: modele, carele
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_ErrorIndic), intent(in) :: ds_errorindic
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_InOut), intent(inout) :: ds_inout
    integer(kind=8) :: fonact(*)
    type(NL_DS_Print), intent(in) :: ds_print
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Init
!
! Prepare storing
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  FONACT : FONCTIONNALITES ACTIVEES
! In  ds_print         : datastructure for printing parameters
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  SDDYNA : SD DYNAMIQUE
! IO  ds_inout         : datastructure for input/output management
! IN  SDCRIT : INFORMATIONS RELATIVES A LA CONVERGENCE
! IN  SDPILO : SD PILOTAGE
! IO  ds_measure       : datastructure for measure and statistics management
! IN  SDERRO : SD ERREUR
! In  ds_energy        : datastructure for energy management
! In  ds_errorindic    : datastructure for error indicator
! IN  MODELE : NOM DU MODELE
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: sdarch
    character(len=24) :: sdarchAinfJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    integer(kind=8) :: numeStoring, numeInst
    integer(kind=8) :: ifm, niv
    aster_logical :: lreuse
    character(len=8) :: result
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_25')
    end if
!
! --- FONCTIONNALITES ACTIVEES
!
    lreuse = isfonc(fonact, 'REUSE')

! - Initial state
    numeInst = 0

! - Get name of result's datastructure
    result = ds_inout%result

! - Name of datastructures
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'

! - Current storing index
    call jeveuo(sdarchAinfJv, 'L', vi=sdarchAinf)
    numeStoring = sdarchAinf(1)

! - Create new datastructure
    if (lreuse) then
        ASSERT(numeStoring .ne. 0)
        call rsrusd(result, numeStoring)
    else
        ASSERT(numeStoring .eq. 0)
        call rscrsd('G', result, 'EVOL_NOLI', 100)
    end if

! - Save initial state
    if (.not. lreuse) then
        call utmess('I', 'ARCHIVAGE_4')
        call nmarch(numeInst, modele, ds_material, carele, fonact, &
                    ds_print, sddisc, sdcrit, &
                    ds_measure, sderro, sddyna, sdpilo, ds_energy, &
                    ds_inout, ds_errorindic, lStoringInitState_=ASTER_TRUE)
    end if
!
end subroutine
