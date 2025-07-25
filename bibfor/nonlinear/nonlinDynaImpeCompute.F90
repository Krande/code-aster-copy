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
! person_in_lload_name: mickael.abbas at edf.fr
!
subroutine nonlinDynaImpeCompute(phaseType, sddyna, &
                                 model, nume_dof, &
                                 ds_material, ds_measure, &
                                 hval_incr, &
                                 hval_veelem, hval_veasse)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asasve.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmtime.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
#include "asterfort/veimpd.h"
!
    integer(kind=8), intent(in) :: phaseType
    character(len=19), intent(in) :: sddyna
    character(len=24), intent(in) :: model, nume_dof
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: hval_veelem(*), hval_veasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Dynamic - Compute impedance
!
! --------------------------------------------------------------------------------------------------
!
! In  phaseType        : name of current phase of algorithm
! In  sddyna           : datastructure for dynamic
! In  model            : name of model
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  ds_material      : datastructure for material parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_veelem      : hat-variable for elementary vectors
! In  hval_veasse      : hat-variable for vectors (node fields)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: vect_elem, vect_asse
    character(len=24) :: vect_alem
    character(len=19) :: vite_prev, vite_curr
    character(len=24), pointer :: v_vect_alem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_6')
    end if
!
! - Get hat variables
!
    call nmchex(hval_incr, 'VALINC', 'VITMOI', vite_prev)
    call nmchex(hval_incr, 'VALINC', 'VITPLU', vite_curr)
    call nmchex(hval_veelem, 'VEELEM', 'CNIMPE', vect_elem)
    call nmchex(hval_veasse, 'VEASSE', 'CNIMPE', vect_asse)
!
! - Launch timer
!
    call nmtime(ds_measure, 'Init', '2nd_Member')
    call nmtime(ds_measure, 'Launch', '2nd_Member')
!
! - Compute
!
    if (phaseType .eq. PRED_EULER) then
        call veimpd(model, ds_material%mateco, vite_prev, sddyna, &
                    vect_elem)
    elseif (phaseType .eq. CORR_NEWTON) then
        call veimpd(model, ds_material%mateco, vite_curr, sddyna, &
                    vect_elem)
    else
        ASSERT(ASTER_FALSE)
    end if
    call asasve(vect_elem, nume_dof, 'R', vect_alem)
    call jeveuo(vect_alem, 'L', vk24=v_vect_alem)
    call copisd('CHAMP_GD', 'V', v_vect_alem(1), vect_asse)
!
! - Stop timer
!
    call nmtime(ds_measure, 'Stop', '2nd_Member')
!
! - Debug
!
    if (niv .ge. 2) then
        call nmdebg('VECT', vect_asse, 6)
    end if
!
end subroutine
