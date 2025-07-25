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
    subroutine vpstur(lmatk, valshi, lmatm, lmatsh, mantis,&
                      expo, pivot, ier, solveu, caldet,&
                      calfac)
        integer(kind=8) :: lmatk
        real(kind=8) :: valshi
        integer(kind=8) :: lmatm
        integer(kind=8) :: lmatsh
        real(kind=8) :: mantis
        integer(kind=8) :: expo
        integer(kind=8) :: pivot
        integer(kind=8) :: ier
        character(len=19) :: solveu
        aster_logical :: caldet
        aster_logical :: calfac
    end subroutine vpstur
end interface
