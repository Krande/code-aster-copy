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
subroutine nonlinDynaMDampCompute(phaseType, &
                                  nlDynaDamping, &
                                  nume_dof, ds_measure, &
                                  hval_incr, hval_veasse)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/fmodam.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmtime.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: phaseType
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: hval_veasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Dynamic - Compute modal damping
!
! --------------------------------------------------------------------------------------------------
!
! In  phaseType        : name of current phase of algorithm
! In  nlDynaDamping    : damping parameters
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
    character(len=19) :: vectAsse, viteIter, viteCurr
    real(kind=8), pointer :: valeVectAsse(:) => null()
    integer(kind=8) :: nb_equa
    type(MODAL_DAMPING) :: modalDamping
    character(len=24) :: jvDataDamp
    character(len=24) :: valmod, basmod
    aster_logical :: lReacVite
    real(kind=8), pointer :: valeViteIter(:) => null()
    real(kind=8), pointer :: valeViteCurr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_10')
    end if

! - Get hat variables
    call nmchex(hval_incr, 'VALINC', 'VITPLU', viteCurr)
    call nmchex(hval_incr, 'VALINC', 'VITKM1', viteIter)
    call nmchex(hval_veasse, 'VEASSE', 'CNAMOD', vectAsse)
    call jeveuo(viteIter(1:19)//'.VALE', 'E', vr=valeViteIter)
    call jeveuo(viteCurr(1:19)//'.VALE', 'E', vr=valeViteCurr)
    call jeveuo(vectAsse(1:19)//'.VALE', 'E', vr=valeVectAsse)

! - Initializations
    call dismoi('NB_EQUA', nume_dof, 'NUME_DDL', repi=nb_equa)
    modalDamping = nlDynaDamping%modalDamping
    jvDataDamp = modalDamping%jvDataDamp
    lReacVite = modalDamping%lReacVite
    valmod = jvDataDamp(1:19)//'.VALM'
    basmod = jvDataDamp(1:19)//'.BASM'

! - Launch timer
    call nmtime(ds_measure, 'Init', '2nd_Member')
    call nmtime(ds_measure, 'Launch', '2nd_Member')

! - Compute
    if (phaseType .eq. PRED_EULER) then
        call fmodam(nb_equa, valeViteIter, valmod, basmod, valeVectAsse)
    elseif (phaseType .eq. CORR_NEWTON) then
        call fmodam(nb_equa, valeViteCurr, valmod, basmod, valeVectAsse)
        if (lReacVite) then
            call fmodam(nb_equa, valeViteCurr, valmod, basmod, valeVectAsse)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if

! - Stop timer
    call nmtime(ds_measure, 'Stop', '2nd_Member')

! - Debug
    if (niv .ge. 2) then
        call nmdebg('VECT', vectAsse, 6)
    end if
!
end subroutine
