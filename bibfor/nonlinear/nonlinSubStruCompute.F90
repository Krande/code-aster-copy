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
subroutine nonlinSubStruCompute(ds_measure, disp, &
                                hval_measse, cnsstr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmtime.h"
#include "asterfort/nmmacv.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/nmdebg.h"
!
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: disp
    character(len=19), intent(in) :: cnsstr, hval_measse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute sub-structuring effect on second member
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_measure       : datastructure for measure and statistics management
! In  disp             : displacement field
! In  hval_measse      : hat-variable for matrix
! In  cnsstr           : name of sub-structuring effect on second member
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: matr_sstr
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_7')
    end if
!
! - Get hat variables
!
    call nmchex(hval_measse, 'MEASSE', 'MESSTR', matr_sstr)
!
! - Launch timer
!
    call nmtime(ds_measure, 'Init', '2nd_Member')
    call nmtime(ds_measure, 'Launch', '2nd_Member')
!
! - Compute
!
    call nmmacv(disp, matr_sstr, cnsstr)
!
! - Stop timer
!
    call nmtime(ds_measure, 'Stop', '2nd_Member')
!
! - Debug
!
    if (niv .ge. 2) then
        call nmdebg('VECT', cnsstr, 6)
    end if
!
end subroutine
