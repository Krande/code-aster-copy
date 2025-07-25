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
subroutine dbrReadOrtho(paraOrtho)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_ParaDBR_Ortho), intent(inout) :: paraOrtho
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Read parameters - For orthogonalization
!
! --------------------------------------------------------------------------------------------------
!
! IO  paraOrtho        : datastructure for parameters (orthogonalization)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nocc
    character(len=8) :: baseInitName
    real(kind=8) :: alpha
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_4')
    end if
!
! - Initializations
!
    baseInitName = ' '
    alpha = 0.d0
!
! - Get parameters
!
    call getvid(' ', 'BASE', scal=baseInitName, nbret=nocc)
    if (nocc .eq. 0) then
        baseInitName = ' '
    end if
    call getvr8(' ', 'ALPHA', scal=alpha, nbret=nocc)
    ASSERT(nocc .eq. 1)
!
! - Save parameters in datastructure
!
    paraOrtho%alpha = alpha
    paraOrtho%baseInitName = baseInitName
!
end subroutine
