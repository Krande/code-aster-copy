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
subroutine dbrChckOrtho(paraOrtho, lReuse)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/romModeChck.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_ParaDBR_Ortho), intent(in) :: paraOrtho
    aster_logical, intent(in) :: lReuse
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Some checks - For orthogonalization
!
! --------------------------------------------------------------------------------------------------
!
! In  paraOrtho        : datastructure for parameters (orthogonalization)
! In  lReuse           : .true. if reuse
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    type(ROM_DS_Field) :: mode
    character(len=8) :: baseInitName
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_35')
    end if
!
! - Initialisations
!
    mode = paraOrtho%baseInit%mode
!
! - Check empiric mode
!
    if (.not. lReuse) then
        call romModeChck(mode)
    end if
!
! - No reuse:
!
    baseInitName = paraOrtho%baseInitName
    if (lReuse) then
        if (baseInitName .ne. ' ') then
            call utmess('F', 'ROM18_21')
        end if
    end if
!
! - Only on nodal fields
!
    if (mode%fieldSupp .ne. 'NOEU') then
        call utmess('F', 'ROM18_36')
    end if
!
end subroutine
