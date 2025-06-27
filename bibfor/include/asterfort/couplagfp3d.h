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
#include "asterf_types.h"
interface 
       subroutine couplagfp3d(a, ngf, na, nc,&
                              dpfa_ds, dpfa_dpg, dpg_depsa6, raideur66p)
        real(kind=8), intent(inout) :: a(ngf,ngf+1)
        integer(kind=8), intent(in) :: ngf
        integer(kind=8), intent(in) :: na
        integer(kind=8), intent(in) :: nc
        real(kind=8), intent(in) :: dpfa_ds(nc, 6)
        real(kind=8), intent(in) :: dpfa_dpg(nc)
        real(kind=8), intent(in) :: dpg_depsa6(6)
        real(kind=8), intent(in) :: raideur66p(6, 6)
    end subroutine couplagfp3d
end interface 
