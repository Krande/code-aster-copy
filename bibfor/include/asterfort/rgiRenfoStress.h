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
    subroutine rgiRenfoStress(xmat, iadrmat, sigmf6, epstf6, epspt6,&
                          teta1, teta2, dt, ppas, theta, fl3d,&
                          end3d, wpl3, vwpl33, vwpl33t, dt3, dr3, ipzero,&
                          ngf, rc00, var0, varf, sigf6d, matdechac, rhov, ierr1)
        integer(kind=8), intent(in) :: ngf
        real(kind=8), intent(in) :: dt3(3)
        real(kind=8), intent(in) :: dr3(3)
        integer(kind=8), intent(in) :: iadrmat
        integer(kind=8), intent(in) :: ipzero(ngf)
        real(kind=8), intent(in) :: xmat(*)
        real(kind=8), intent(in) :: sigmf6(6)
        real(kind=8), intent(in) :: epstf6(6)
        real(kind=8), intent(in) :: epspt6(6)
        real(kind=8), intent(in) :: var0(*)
        real(kind=8), intent(in) :: rc00
        real(kind=8), intent(in) :: teta1
        real(kind=8), intent(in) :: teta2
        real(kind=8), intent(in) :: dt
        real(kind=8), intent(in) :: theta
        real(kind=8), intent(in) :: wpl3(3)
        real(kind=8), intent(in) :: vwpl33(3,3)
        real(kind=8), intent(in) :: vwpl33t(3,3)
        aster_logical, intent(in) :: end3d
        aster_logical, intent(in) :: fl3d
        aster_logical, intent(in) :: ppas
        real(kind=8), intent(out) :: varf(*)
        real(kind=8), intent(out) :: sigf6d(6)
        real(kind=8), intent(out) :: matdechac(6,6)
        real(kind=8), intent(out) :: rhov
        integer(kind=8), intent(out) :: ierr1
    end subroutine rgiRenfoStress
end interface
