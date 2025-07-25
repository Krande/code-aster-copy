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

subroutine matrot(angl_naut, pgl)
!
    implicit none
!
!
    real(kind=8), intent(in) :: angl_naut(*)
    real(kind=8), intent(out) :: pgl(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
! Reference frames management
!
! Construct matrix from nautical angles
!
! --------------------------------------------------------------------------------------------------
!
! From (X0,Y0,Z0) to (X,Y,Z)
!
!       (X0,Y0,Z0)     >    (X1,Y1,Z0)    >    (X,Y1,Z2)    >    (X,Y,Z)
!                    APLHA              BETA              GAMMA
!
! In  angl_naut        : nautical angles
!                        (1) Alpha - clockwise around Z0
!                        (2) Beta  - counterclockwise around Y1
!                        (1) Gamma - clockwise around X
! Out pgl              : matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: cosa, cosb, cosg, sina, sinb, sing
!
! --------------------------------------------------------------------------------------------------
!
    cosa = cos(angl_naut(1))
    sina = sin(angl_naut(1))
    cosb = cos(angl_naut(2))
    sinb = sin(angl_naut(2))
    cosg = cos(angl_naut(3))
    sing = sin(angl_naut(3))
!
    pgl(1, 1) = cosb*cosa
    pgl(2, 1) = sing*sinb*cosa-cosg*sina
    pgl(3, 1) = sing*sina+cosg*sinb*cosa
    pgl(1, 2) = cosb*sina
    pgl(2, 2) = cosg*cosa+sing*sinb*sina
    pgl(3, 2) = cosg*sinb*sina-cosa*sing
    pgl(1, 3) = -sinb
    pgl(2, 3) = sing*cosb
    pgl(3, 3) = cosg*cosb
!
end subroutine
