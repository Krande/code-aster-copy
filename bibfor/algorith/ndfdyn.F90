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
subroutine ndfdyn(sddyna, nlDynaDamping, &
                  hval_incr, hval_measse, ds_measure, &
                  cndyna)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use NonLinearDyna_module, only: compAcceForce, compViteForce
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmtime.h"
#include "asterfort/utmess.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
!
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19), intent(in) :: hval_incr(*), hval_measse(*)
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: cndyna
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! CALCUL DES FORCES DE RAPPEL DYNAMIQUE
!
! --------------------------------------------------------------------------------------------------
!
! In  nlDynaDamping    : damping parameters
! In  sddyna           : datastructure for dynamic
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_measse      : hat-variable for matrix
! IO  ds_measure       : datastructure for measure and statistics management
! In  cndyna           : name of dynamic effect on second member
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19), parameter :: cniner = '&&CNPART.CHP1', cnhyst = '&&CNPART.CHP2'
    real(kind=8) :: coefIner, coefDamp
    aster_logical :: lDampMatrix, l_impl
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_9')
    end if

! - Launch timer
    call nmtime(ds_measure, 'Init', '2nd_Member')
    call nmtime(ds_measure, 'Launch', '2nd_Member')

! - Coefficients
    coefIner = ndynre(sddyna, 'COEF_FDYN_MASSE')
    coefDamp = ndynre(sddyna, 'COEF_FDYN_AMORT')

! - Active functionnalities
    lDampMatrix = nlDynaDamping%hasMatrDamp
    l_impl = ndynlo(sddyna, 'IMPLICITE')

! - Initializations
    call vtzero(cndyna)

! - Compute
    if (l_impl) then
        call compAcceForce(hval_incr, hval_measse, cniner)
        call vtaxpy(coefIner, cniner, cndyna)
    end if
    if (lDampMatrix) then
        call compViteForce(nlDynaDamping, hval_incr, 'VITPLU', cnhyst)
        call vtaxpy(coefDamp, cnhyst, cndyna)
    end if

! - Stop timer
    call nmtime(ds_measure, 'Stop', '2nd_Member')

! - Debug
    if (niv .ge. 2) then
        call nmdebg('VECT', cndyna, 6)
    end if
!
end subroutine
