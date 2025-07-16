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
! person_in_charge: matthieu-m.le-cren at edf.fr
!

subroutine hatSmooth(nno, nnos, v_basf, vector)
!
    implicit none
!
    integer(kind=8), intent(in) :: nno, nnos
    real(kind=8), intent(in), pointer, dimension(:) :: v_basf
    real(kind=8), intent(inout), dimension(nno) :: vector
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    smooth vector G or K in quadratic, using hat shape functions
!
!---------------------------------------------------------------------------------------------------
    !integer :: nnos
    real(kind=8), dimension(nnos-1) :: le
    real(kind=8), dimension(nnos-2) :: lg, ld
    real(kind=8), dimension(nnos) :: smooth
!
    le = abs(v_basf(3:nno:2)-v_basf(1:nno-2:2))
    lg = 2.*le(1:nnos-2)/(le(1:nnos-2)+le(2:nnos-1))
    ld = 2.*le(2:nnos-1)/(le(1:nnos-2)+le(2:nnos-1))
    smooth(1) = (2*vector(1)+vector(2))/3
    smooth(2:nnos-1) = (lg*vector(2:nno-3:2)+vector(3:nno-2:2)+ld*vector(4:nno-1:2))/3
    smooth(nnos) = (vector(nno-1)+2*vector(nno))/3
    vector(1:nno:2) = smooth
    vector(2:nno-1:2) = (smooth(1:nnos-1)+smooth(2:nnos))/2
!
end subroutine
