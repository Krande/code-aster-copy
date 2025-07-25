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
subroutine ncpact(dg, nec, nbc)
    implicit none
!
!
#include "asterfort/nbbit1.h"
    integer(kind=8) :: dg(*), nec, nbc
!
!**********************************************************************
!
!     NBC := NBR DE COMPOSSANTES REPRESENTANT LA GRANDEUR DE
!            DESCRIPTEUR DG ( TAILLE = NEC)
!
!**********************************************************************
!
    integer(kind=8) :: i, compt
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    nbc = 0
    compt = 0
!
    do i = 1, nec, 1
!
        call nbbit1(dg(i), compt)
!
        nbc = nbc+compt
!
    end do
!
end subroutine
