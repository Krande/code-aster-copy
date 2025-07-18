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
    subroutine mmproj(alias, nno, ndim, coorma, coorpt,&
                      itemax, epsmax, toleou, dirapp, dir,&
                      ksi1, ksi2, tau1, tau2, iproj,&
                      niverr)
        character(len=8) :: alias
        integer(kind=8) :: nno
        integer(kind=8) :: ndim
        real(kind=8) :: coorma(27)
        real(kind=8) :: coorpt(3)
        integer(kind=8) :: itemax
        real(kind=8) :: epsmax
        real(kind=8) :: toleou
        aster_logical :: dirapp
        real(kind=8) :: dir(3)
        real(kind=8) :: ksi1
        real(kind=8) :: ksi2
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
        integer(kind=8) :: iproj
        integer(kind=8) :: niverr
    end subroutine mmproj
end interface
