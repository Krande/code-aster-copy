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

subroutine getElemOrientation(ndim, nno, jv_geom, angl_naut)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/rcangm.h"
!
    integer(kind=8), intent(in) :: ndim, nno
    integer(kind=8), intent(in) :: jv_geom
    real(kind=8), intent(out) :: angl_naut(3)
!
! --------------------------------------------------------------------------------------------------
!
!
! Get element frame orientation for anisotropy
!
! --------------------------------------------------------------------------------------------------
!
! In  ndim             : dimension of element (2 ou 3)
! In  nno              : number of nodes (all)
! In  jv_geom          : JEVEUX adress to initial geometry (mesh)
! Out angl_naut        : nautical angles for frame orientation
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: coor(3)
    integer(kind=8) :: i_node, i_dim
!
! --------------------------------------------------------------------------------------------------
!
    coor(:) = 0.d0
    angl_naut(:) = 0.d0
!
! - Compute barycentric center
!
    do i_node = 1, nno
        do i_dim = 1, ndim
            coor(i_dim) = coor(i_dim)+zr(jv_geom+i_dim+ndim*(i_node-1)-1)/nno
        end do
    end do
!
! - Get nautical angles
!
    call rcangm(ndim, coor, angl_naut)
!
end subroutine
