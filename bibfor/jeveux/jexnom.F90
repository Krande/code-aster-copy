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

function jexnom(nomc, nomo)
    implicit none
#include "jeveux.h"
!
    character(len=32) :: jexnom
    character(len=*), intent(in) :: nomc, nomo
!     ------------------------------------------------------------------
    integer(kind=8) :: numec
    common/inumje/numec
    real(kind=8) :: reelc
    common/reelje/reelc
    character(len=24) :: nomec
    common/knomje/nomec
!     ------------------------------------------------------------------
    character(len=24) :: ch24
    character(len=8) :: ch8
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    data ch8/'$$XNOM  '/
!     ------------------------------------------------------------------
!
    numec = 0
    reelc = 0.d0
    nomec = nomo
    ch24 = nomc
    jexnom(1:24) = ch24
    jexnom(25:32) = ch8
end function
