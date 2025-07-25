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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine mmgeom(ndim, &
                  nne, nnm, &
                  ffe, ffm, &
                  elem_slav_coor, elem_mast_coor, &
                  tau1, tau2, &
                  norm, mprojn, mprojt, &
                  geome, geomm)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/mmcaln.h"
!
    integer(kind=8), intent(in) :: ndim, nne, nnm
    real(kind=8), intent(in) :: ffe(9), ffm(9)
    real(kind=8), intent(in) :: elem_slav_coor(9, 3), elem_mast_coor(9, 3)
    real(kind=8), intent(in) :: tau1(3), tau2(3)
    real(kind=8), intent(out) :: norm(3), mprojn(3, 3), mprojt(3, 3)
    real(kind=8), intent(out) :: geomm(3), geome(3)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Compute local basis
!
! --------------------------------------------------------------------------------------------------
!
! In  ndim             : dimension of problem (2 or 3)
! In  nne              : number of slave nodes
! In  nnm              : number of master nodes
! In  ffe              : shape function for slave nodes
! In  ffm              : shape function for master nodes
! In  elem_slav_coor   : updated coordinates from slave side of contact element
! In  elem_mast_coor   : updated coordinates from master side of contact element
! In  tau1             : first tangent at current contact point
! In  tau2             : second tangent at current contact point
! Out norm             : normal at current contact point
! Out mprojn           : matrix of normal projection
! Out mprojt           : matrix of tangent projection
! Out geome            : coordinates for contact point
! Out geomm            : coordinates for projection of contact point
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: idim, inoe, inom
!
! --------------------------------------------------------------------------------------------------
!
    geome(:) = 0.d0
    geomm(:) = 0.d0
!
! - Coordinates for contact point
!
    do idim = 1, ndim
        do inoe = 1, nne
            geome(idim) = geome(idim)+ffe(inoe)*elem_slav_coor(inoe, idim)
        end do
    end do
!
! - Coordinates for projection of contact point
!
    do idim = 1, ndim
        do inom = 1, nnm
            geomm(idim) = geomm(idim)+ffm(inom)*elem_mast_coor(inom, idim)
        end do
    end do
!
! - Compute local basis
!
    call mmcaln(ndim, tau1, tau2, norm, mprojn, mprojt)
!
end subroutine
