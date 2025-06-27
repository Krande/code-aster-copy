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
    subroutine rc36f1(nbsigr, nocc, saltij, isk, isl,&
                      nk, nl, n0, nbp12, nbp23,&
                      nbp13, sigr, yapass, typass, nsitup)
        integer(kind=8) :: nbsigr
        integer(kind=8) :: nocc(*)
        real(kind=8) :: saltij(*)
        integer(kind=8) :: isk
        integer(kind=8) :: isl
        integer(kind=8) :: nk
        integer(kind=8) :: nl
        integer(kind=8) :: n0
        integer(kind=8) :: nbp12
        integer(kind=8) :: nbp23
        integer(kind=8) :: nbp13
        integer(kind=8) :: sigr(*)
        aster_logical :: yapass
        character(len=3) :: typass
        integer(kind=8) :: nsitup
    end subroutine rc36f1
end interface
