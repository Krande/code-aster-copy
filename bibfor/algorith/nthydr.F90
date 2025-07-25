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
subroutine nthydr(hydrat)
    implicit none
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/getvtx.h"
    aster_logical :: hydrat
    integer(kind=8) :: nbocc, n1, i
    character(len=16) :: comp
!     ------------------------------------------------------------------
!
    hydrat = .false.
!
    call getfac('COMPORTEMENT', nbocc)
!
    do i = 1, nbocc
!
        call getvtx('COMPORTEMENT', 'RELATION', iocc=i, scal=comp, nbret=n1)
!
        if (comp(1:9) .eq. 'THER_HYDR') hydrat = .true.
!
    end do
!
! FIN ------------------------------------------------------------------
end subroutine
