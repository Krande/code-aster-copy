! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine matrHookePlaneStrain(elas_type, repere,&
                                h, g, g1, &
                                matr_elas)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dpao2d.h"
#include "asterfort/utbtab.h"
!
!
    integer, intent(in) :: elas_type
    real(kind=8), intent(in) :: repere(7)
    real(kind=8), intent(in) :: h(6), g
    real(kind=8), intent(in) :: g1
    real(kind=8), intent(out) :: matr_elas(4, 4)
!
! --------------------------------------------------------------------------------------------------
!
! Hooke matrix for iso-parametric elements
!
! Plane strain and axisymmetric
!
! --------------------------------------------------------------------------------------------------
!
! In  elas_id          : Type of elasticity
!                 1 - Isotropic
!                 2 - Orthotropic
!                 3 - Transverse isotropic
!                           or viscoelasticity
!                 4 - Isotropic
!                 5 - Orthotropic
!                 6 - Transverse isotropic
! In  repere           : define reference frame (AFFE_CARA_ELEM/MASSIF)
!                        repere(1) =  1 => nautical angles (ANGL_REP)
!                           repere(2:4) : nautical angles
!                           repere(5:7) : 0.d0
!                        repere(1) =  2 => Euler angles (ANGL_EULER)
!                           repere(2:4) : nautical angles
!                           repere(5:7) : Euler angles
!                        repere(1) = -1 => axisymetric axis (ANGL_AXE)
!                           repere(2:4) : ANGL_AXE
!                           repere(5:7) : ORIG_AXE
! In  h                : Hook coefficient (all)
! In  g                : shear ratio (isotropic/Transverse isotropic)
! In  g1               : shear ratio (Orthotropic)
! Out matr_elas        : Hooke matrix
! In  xyzgau           : coordinates of Gauss point (for ANGL_AXE case)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: irep, i, j
    real(kind=8) :: matr_tran(4, 4), dorth(4, 4), work(4, 4)
!
! --------------------------------------------------------------------------------------------------
!
    matr_elas(:,:) = 0.d0
    dorth(:,:)     = 0.d0
    work(:,:)      = 0.d0
!
! - Compute Hooke matrix
!
    if (elas_type .eq. 1 .or. elas_type .eq. 4) then
!
! ----- Isotropic matrix
!
        matr_elas(1,1) = h(1)
        matr_elas(1,2) = h(2)
        matr_elas(1,3) = h(2)
        matr_elas(2,1) = h(2)
        matr_elas(2,2) = h(1)
        matr_elas(2,3) = h(2)
        matr_elas(3,1) = h(2)
        matr_elas(3,2) = h(2)
        matr_elas(3,3) = h(1)
        matr_elas(4,4) = g
!
    else if (elas_type .eq. 2 .or. elas_type .eq. 5) then
!
! ----- Orthotropic matrix
!
        dorth(1,1) = h(1)
        dorth(1,2) = h(2)
        dorth(1,3) = h(3)
        dorth(2,2) = h(4)
        dorth(2,3) = h(5)
        dorth(3,3) = h(6)
        dorth(2,1) = dorth(1,2)
        dorth(3,1) = dorth(1,3)
        dorth(3,2) = dorth(2,3)
        dorth(4,4) = g1
!
! ----- Matrix from orthotropic basis to global 3D basis
!
        call dpao2d(repere, irep, matr_tran)
!
! ----- Hooke matrix in global 3D basis
!
        ASSERT((irep.eq.1).or.(irep.eq.0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 4, 4, dorth, matr_tran, work, matr_elas)
        else if (irep.eq.0) then
            do i = 1, 4
                do j = 1, 4
                    matr_elas(i,j) = dorth(i,j)
                end do
            end do
        endif
!
    else if (elas_type .eq. 3 .or. elas_type .eq. 6) then
!
! ----- Transverse isotropic matrix
!
        dorth(1,1) = h(1)
        dorth(1,2) = h(2)
        dorth(1,3) = h(3)
        dorth(2,1) = dorth(1,2)
        dorth(2,2) = dorth(1,1)
        dorth(2,3) = dorth(1,3)
        dorth(3,1) = dorth(1,3)
        dorth(3,2) = dorth(2,3)
        dorth(3,3) = h(4)
        dorth(4,4) = h(5)
!
! ----- Matrix from transverse isotropic basis to global 3D basis
!
        call dpao2d(repere, irep, matr_tran)
!
! ----- Hooke matrix in global 3D basis
!
        ASSERT((irep.eq.1).or.(irep.eq.0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 4, 4, dorth, matr_tran, work, matr_elas)
        else if (irep.eq.0) then
            do i = 1, 4
                do j = 1, 4
                    matr_elas(i,j) = dorth(i,j)
                end do
            end do
        endif
    else
        ASSERT(.false.)
    endif
!
end subroutine
