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
    subroutine apm345(nbtetc, typcon, rayonc, centrc, nk,&
                      k24rc, pivot2, ltest, typcha, lraide,&
                      lmasse, ldynam, solveu, lamor, lc,&
                      impr, ifapm)
        integer(kind=8) :: nbtetc
        character(len=8) :: typcon
        real(kind=8) :: rayonc
        complex(kind=8) :: centrc
        integer(kind=8) :: nk
        character(len=24) :: k24rc
        integer(kind=8) :: pivot2
        aster_logical :: ltest
        character(len=8) :: typcha
        integer(kind=8) :: lraide
        integer(kind=8) :: lmasse
        integer(kind=8) :: ldynam
        character(len=19) :: solveu
        integer(kind=8) :: lamor
        aster_logical :: lc
        character(len=3) :: impr
        integer(kind=8) :: ifapm
    end subroutine apm345
end interface
