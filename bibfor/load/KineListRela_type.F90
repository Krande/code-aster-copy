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
! ==================================================================================================
!
! Types to manage list of kinematic relations
!
! ==================================================================================================
!
module KineListRela_type
!
    implicit none
! ==================================================================================================
#include "asterf_types.h"
! ==================================================================================================
    type KINE_LIST_RELA
! - Type of relation: (implicit or explicit)
!       'implicit': no second member (RHS)
!       'explicit': with a second membre (RHS)
        character(len=16) :: relaType = ' '

! - Number of terms in relation
        integer(kind=8) :: nbTermMaxi = 0

! - List of relations (LHS)
        character(len=4) :: coefMultType = ' '
        character(len=8), pointer :: nodeName(:) => null()
        character(len=8), pointer :: dofName(:) => null()
        real(kind=8), pointer :: coefMultReal(:) => null()
        real(kind=8) :: coefMultTole = 1.d-6

! - List of relations (RHS)
        character(len=4) :: coefImpoType = ' '
        real(kind=8) :: coefImpoReal = 0.d0
        character(len=8) :: coefImpoFunc = ' '
        complex(kind=8) :: coefImpoCplx = dcmplx(0.d0, 0.d0)

! - Local coordinates system
        real(kind=8), pointer :: LCSVale(:) => null()
        integer(kind=8), pointer :: LCSType(:) => null()

! - Name of JEVEUX object for aflrch subroutine
        character(len=19) :: listLineRela = ' '

    end type KINE_LIST_RELA
!
end module KineListRela_type
