! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
module HHO_inertia_module
!
    use HHO_type
    use HHO_quadrature_module
    use HHO_geometry_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "blas/dsyev.h"
#include "blas/dsyr.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - utilitaries
!
! Some common utilitaries
!
! --------------------------------------------------------------------------------------------------
    public :: hhoLocalAxesCell, hhoLocalAxesFace
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoLocalAxesCell(hhoCell) result(axes)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8) :: axes(3, 3)
!
! ---------------------------------------------------------------------------------
!  HHO - inertial
!  Compute inertia axes
!
!   In hhoCell     : the current HHO Cell
! ---------------------------------------------------------------------------------
!
        type(HHO_quadrature)  :: hhoQuad
        integer :: ipg, idim
        integer(kind=4) :: info
        real(kind=8) :: coor(3), evalues(3), work(50)
!
        axes = 0.d0
!
!
! ----- get quadrature
        call hhoQuad%GetQuadCell(hhoCell, 2)
!
! ----- Loop on quadrature point
        do ipg = 1, hhoQuad%nbQuadPoints
            coor = hhoCell%barycenter-hhoQuad%points(1:3, ipg)
            call dsyr('U', hhoCell%ndim, hhoQuad%weights(ipg), coor, 1, axes, 3)
        end do
!
! ----- Compute eigenvector
        evalues = 0.d0
        call dsyev('V', 'U', hhoCell%ndim, axes, 3, evalues, work, 50, info)
        ASSERT(info == 0)
!
        do idim = 1, hhoCell%ndim
            axes(1:3, idim) = axes(1:3, idim)/norm2(axes(1:3, idim))
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoLocalAxesFace(hhoFace) result(axes)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        real(kind=8) :: axes(3, 2)
!
! ---------------------------------------------------------------------------------
!  HHO - inertial
!  Compute inertia axes
!
!   In hhoFace     : the current HHO Face
! ---------------------------------------------------------------------------------
!
        type(HHO_quadrature)  :: hhoQuad
        integer :: ipg, idim
        integer(kind=4) :: info
        real(kind=8) :: coor(3), evalues(3), work(50), axes_3d(3, 3)
!
        axes = 0.d0
        axes_3d = 0.d0
!
        if (hhoFace%ndim == 1) then
            coor = hhoFace%coorno(1:3, 2)-hhoFace%coorno(1:3, 1)
            axes(1:3, 1) = coor/norm2(coor)
        else
!
! ----- get quadrature
            call hhoQuad%GetQuadFace(hhoFace, 2)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
                coor = hhoFace%barycenter-hhoQuad%points(1:3, ipg)
                call dsyr('U', hhoFace%ndim+1, hhoQuad%weights(ipg), coor, 1, axes_3d, 3)
            end do
!
! ----- Compute eigenvector
            evalues = 0.d0
            call dsyev('V', 'U', hhoFace%ndim+1, axes_3d, 3, evalues, work, 50, info)
            ASSERT(info == 0)
            ASSERT(minloc(evalues(1:hhoFace%ndim+1), dim=1) == 1)
            axes(1:3, 1:2) = axes_3d(1:3, 2:3)
!
            do idim = 1, hhoFace%ndim
                axes(1:3, idim) = axes(1:3, idim)/norm2(axes(1:3, idim))
            end do
        end if
!
    end function
!
end module
