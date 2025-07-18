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
subroutine nmini0(eta, nume_inst, matass, &
                  zmeelm, zmeass, zveelm, &
                  zveass, zsolal, zvalin, &
                  ds_print, ds_conv, ds_algopara, &
                  ds_inout, ds_contact, ds_measure, &
                  ds_energy, ds_material, sderro)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/nmchai.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/nonlinDSConvergenceCreate.h"
#include "asterfort/nonlinDSPrintCreate.h"
#include "asterfort/nonlinDSAlgoParaCreate.h"
#include "asterfort/nonlinDSInOutCreate.h"
#include "asterfort/nonlinDSContactCreate.h"
#include "asterfort/nonlinDSMeasureCreate.h"
#include "asterfort/nonlinDSEnergyCreate.h"
#include "asterfort/nonlinDSMaterialCreate.h"
#include "asterfort/nmcrga.h"
#include "asterfort/nonlinDSPrintSepLine.h"
!
    character(len=19), intent(out) :: matass
    integer(kind=8), intent(out) :: nume_inst
    real(kind=8), intent(out) :: eta
    integer(kind=8), intent(in) :: zmeelm, zmeass, zveelm
    integer(kind=8), intent(in) :: zveass, zsolal, zvalin
    type(NL_DS_Print), intent(out) :: ds_print
    type(NL_DS_Conv), intent(out) :: ds_conv
    type(NL_DS_AlgoPara), intent(out) :: ds_algopara
    type(NL_DS_InOut), intent(out) :: ds_inout
    type(NL_DS_Contact), intent(out) :: ds_contact
    type(NL_DS_Measure), intent(out) :: ds_measure
    type(NL_DS_Energy), intent(out) :: ds_energy
    type(NL_DS_Material), intent(out) :: ds_material
    character(len=24) :: sderro
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Creation of datastructures
!
! --------------------------------------------------------------------------------------------------
!
! Out nume_inst        : index of current time step
! Out ds_print         : datastructure for printing parameters
! Out ds_conv          : datastructure for convergence management
! Out ds_algopara      : datastructure for algorithm parameters
! Out ds_inout         : datastructure for input/output management
! Out ds_contact       : datastructure for contact management
! Out ds_measure       : datastructure for measure and statistics management
! Out ds_energy        : datastructure for energy management
! Out ds_material      : datastructure for material parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: long
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call nonlinDSPrintSepLine()
        call utmess('I', 'MECANONLINE14_99')
    end if
!
! - Create printing management datastructure
!
    call nonlinDSPrintCreate('MECA', ds_print)
!
! - Create convergence management datastructure
!
    call nonlinDSConvergenceCreate(ds_conv)
!
! - Create algorithm parameters datastructure
!
    call nonlinDSAlgoParaCreate(ds_algopara)
!
! - Create input/output management datastructure
!
    call nonlinDSInOutCreate('MECA', ds_inout)
!
! - Create contact management datastructure
!
    call nonlinDSContactCreate(ds_contact)
!
! - Create measure and statistics management datastructure
!
    call nonlinDSMeasureCreate(ds_measure)
!
! - Create energy management datastructure
!
    call nonlinDSEnergyCreate(ds_energy)
!
! - Create material management datastructure
!
    call nonlinDSMaterialCreate(ds_material)

! - Create datastructure for events in algorithm
    call nmcrga(sderro)
!
! --- INITIALISATION BOUCLE EN TEMPS
!
    nume_inst = 0
    eta = zero
    matass = '&&OP0070.MATASS'
!
! --- VERIF. LONGUEURS VARIABLES CHAPEAUX (SYNCHRO OP0070/NMCHAI)
!
    call nmchai('MEELEM', 'LONMAX', long)
    ASSERT(long .eq. zmeelm)
    call nmchai('MEASSE', 'LONMAX', long)
    ASSERT(long .eq. zmeass)
    call nmchai('VEELEM', 'LONMAX', long)
    ASSERT(long .eq. zveelm)
    call nmchai('VEASSE', 'LONMAX', long)
    ASSERT(long .eq. zveass)
    call nmchai('SOLALG', 'LONMAX', long)
    ASSERT(long .eq. zsolal)
    call nmchai('VALINC', 'LONMAX', long)
    ASSERT(long .eq. zvalin)
!
end subroutine
