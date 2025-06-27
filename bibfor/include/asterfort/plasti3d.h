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
    subroutine plasti3d(xmat, inputR, inputVR6, inputMat33, inputI,&
                    inputL, var0, raideur66, souplesse66,&
                    A, B, X, ngf, varf, ipzero,&
                    outputR, outputVR6, outputMat33,&
                    outputI)
        integer(kind=8), intent(in) :: inputI(*)
        integer(kind=8), intent(in) :: ngf
        aster_logical, intent(in) :: inputL(*)
        real(kind=8), intent(in) :: xmat(*)
        real(kind=8), intent(in) :: var0(*)
        real(kind=8), intent(in) :: inputR(*)
        real(kind=8), intent(in) :: inputVR6(6,*)
        real(kind=8), intent(in) :: inputMat33(3,3,*)
        real(kind=8), intent(in) :: raideur66(6,6)
        real(kind=8), intent(in) :: souplesse66(6,6)
        real(kind=8), intent(inout) :: A(ngf, ngf+1)
        real(kind=8), intent(inout) :: B(ngf)
        real(kind=8), intent(inout) :: X(ngf)
        real(kind=8), intent(inout) :: varf(*)
        real(kind=8), intent(out) :: outputR(*)
        real(kind=8), intent(out) :: outputVR6(6,*)
        real(kind=8), intent(out) :: outputMat33(3,3,*)
        integer(kind=8), intent(out) :: outputI(*)
        integer(kind=8), intent(inout) :: ipzero(ngf)
    end subroutine plasti3d
end interface
