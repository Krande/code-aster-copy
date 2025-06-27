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
    subroutine ntweib(nrupt, cals, sk, sigw, nur,&
                      nt, nbres, x1, x2, xacc,&
                      rtsafe, impr, ifm, indtp, nbtp)
        integer(kind=8) :: nrupt
        aster_logical :: cals
        real(kind=8) :: sk(*)
        real(kind=8) :: sigw(*)
        integer(kind=8) :: nur(*)
        integer(kind=8) :: nt(*)
        integer(kind=8) :: nbres
        real(kind=8) :: x1
        real(kind=8) :: x2
        real(kind=8) :: xacc
        real(kind=8) :: rtsafe
        aster_logical :: impr
        integer(kind=8) :: ifm
        integer(kind=8) :: indtp(*)
        integer(kind=8) :: nbtp
    end subroutine ntweib
end interface
