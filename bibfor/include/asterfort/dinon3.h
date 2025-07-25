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
    subroutine dinon3(neq, ul, dul, utl, nno,&
                      nbcomp, varimo, raide, nbpar, param,&
                      okdire, varipl)
        integer(kind=8) :: nbpar
        integer(kind=8) :: nbcomp
        integer(kind=8) :: neq
        real(kind=8) :: ul(neq)
        real(kind=8) :: dul(neq)
        real(kind=8) :: utl(neq)
        integer(kind=8) :: nno
        real(kind=8) :: varimo(nbcomp*3)
        real(kind=8) :: raide(nbcomp)
        real(kind=8) :: param(6, nbpar)
        aster_logical :: okdire(6)
        real(kind=8) :: varipl(nbcomp*3)
    end subroutine dinon3
end interface
