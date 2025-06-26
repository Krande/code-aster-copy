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

subroutine mm_cycl_laugf(pres, dist, coef_augm, lagr_norm)
!
    implicit none
!
! person_in_charge: mickael.abbas at edf.fr
!
    real(kind=8), intent(in) :: pres(3)
    real(kind=8), intent(in) :: dist(3)
    real(kind=8), intent(in) :: coef_augm
    real(kind=8), intent(out) :: lagr_norm
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve - Cycling
!
! Augmented lagrangian (vectorial version)
!
! --------------------------------------------------------------------------------------------------
!
! In  pres      : pressure
! In  dist      : distance
! In  coef_augm : augmented coefficient
! In  ndim      : topological dimension
! Out lagr_norm : norm of augmented lagrangian
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: lagr_augm(3)
!
! --------------------------------------------------------------------------------------------------
!
    lagr_norm = 0.d0
!
! -- Test to prevent FPE
!
    if (maxval(abs(pres)) > 10.d50) then
        lagr_norm = 10.d50
    elseif (maxval(abs(dist)) > 10.d10) then
        lagr_norm = 10.d50
    else
        lagr_augm(1:3) = pres(1:3)+coef_augm*dist(1:3)
        lagr_norm = norm2(lagr_augm)
    end if
!
end subroutine
