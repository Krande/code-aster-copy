! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
    subroutine plastMult(na, nf0, ngf, x, irr, nc, indic2, ipla2,&
                         imax, ig, a, b, dgfa_ds, goto20)
        integer, intent(inout) :: na
        integer, intent(inout) :: ipla2
        integer, intent(inout) :: ig(nc)
        integer, intent(in) :: nf0
        integer, intent(in) :: irr
        integer, intent(in) :: nc
        integer, intent(in) :: imax
        integer, intent(in) :: ngf
        real(kind=8), intent(inout) :: x(ngf)
        real(kind=8), intent(inout) :: a(ngf,ngf+1)
        real(kind=8), intent(inout) :: b(ngf)
        real(kind=8), intent(inout) :: dgfa_ds(nc, 6)
        aster_logical, intent(inout) :: indic2
        aster_logical, intent(out) :: goto20
    end subroutine plastMult
end interface
