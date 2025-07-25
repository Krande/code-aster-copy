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
    subroutine ejcine_hm(ndim, axi, nno1, nno2, vff1,&
                      vff2, wref, dffr2, geom, ang,&
                      wg, b)
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        aster_logical :: axi
        real(kind=8) :: vff1(nno1)
        real(kind=8) :: vff2(nno2)
        real(kind=8) :: wref
        real(kind=8) :: dffr2(ndim-1, nno2)
        real(kind=8) :: geom(ndim, nno2)
        real(kind=8) :: ang(*)
        real(kind=8) :: wg
        real(kind=8) :: b(3, 3, 2*nno1)
    end subroutine ejcine_hm
end interface
