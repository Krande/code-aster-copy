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
subroutine romFSINumberingInit(field, algoGreedy)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/indik8.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_Field), intent(in) :: field
    type(ROM_DS_AlgoGreedy), intent(inout) :: algoGreedy
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Create numbering of nodes for FSI
!
! --------------------------------------------------------------------------------------------------
!
! In  field            : field to analyze
! IO  algoGreedy       : datastructure for Greedy algorithm
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: cmpNume
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_53')
    end if
!
    cmpNume = indik8(field%listCmpName, 'PRES', 1, field%nbCmpName)
    algoGreedy%nume_pres = cmpNume
    cmpNume = indik8(field%listCmpName, 'PHI', 1, field%nbCmpName)
    algoGreedy%nume_phi = cmpNume
!
end subroutine
