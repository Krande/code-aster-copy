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

!
!
#include "asterf_types.h"
!
interface
    subroutine cfnors(noma, ds_contact, posmai, typent,&
                      nument, lpoutr, lpoint, ksi1, ksi2,&
                      lliss, itype, vector, tau1, tau2,&
                      lnfixe)
        use NonLin_Datastructure_type
        character(len=8) :: noma
        type(NL_DS_Contact), intent(in) :: ds_contact
        integer(kind=8) :: posmai
        character(len=4) :: typent
        integer(kind=8) :: nument
        aster_logical :: lpoutr
        aster_logical :: lpoint
        real(kind=8) :: ksi1
        real(kind=8) :: ksi2
        aster_logical :: lliss
        integer(kind=8) :: itype
        real(kind=8) :: vector(3)
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
        aster_logical :: lnfixe
    end subroutine cfnors
end interface
