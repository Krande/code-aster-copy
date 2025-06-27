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
    subroutine tirViscoElas(fl3d, var0, xmat, inputR, inputVR6, ngf,  &
                        deltam, avean, A, B, X, ipzero, &
                        epsk16, epsm16, epse16, &
                        sig16, sigke16, raideur66, we1)
        aster_logical, intent(in) :: fl3d
        real(kind=8), intent(in) :: inputR(*)
        real(kind=8), intent(in) :: inputVR6(6,*)
        real(kind=8), intent(in) :: raideur66(6,6)
        integer(kind=8), intent(in) :: ngf
        integer(kind=8), intent(inout) :: ipzero(ngf)
        real(kind=8), intent(in) :: var0(*)
        real(kind=8), intent(in) :: xmat(*)
        real(kind=8), intent(out) :: deltam
        real(kind=8), intent(out) :: avean
        real(kind=8), intent(out) :: epsk16(6)
        real(kind=8), intent(out) :: epsm16(6)
        real(kind=8), intent(out) :: epse16(6)
        real(kind=8), intent(out) :: sig16(6)
        real(kind=8), intent(out) :: sigke16(6)
        real(kind=8), intent(out) :: we1
        real(kind=8), intent(inout) :: A(ngf, ngf+1)
        real(kind=8), intent(inout) :: B(ngf)
        real(kind=8), intent(inout) :: X(ngf)
    end subroutine tirViscoElas
end interface
