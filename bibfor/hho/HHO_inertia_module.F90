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
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterf_debug.h"
#include "blas/dsyev.h"
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
        type(HHO_quadrature) :: hhoQuad
        integer(kind=8) :: ipg, idim
        real(kind=8) :: coor(3), evalues(3), work(50), qp_coor(3)
        blas_int :: b_lda, b_lwork, b_n, info
!
        axes = 0.d0
!
!
! ----- get quadrature
        ! need approximate quadrature
        call hhoQuad%GetQuadCell(hhoCell, 2, split=ASTER_FALSE)
!
! ----- Loop on quadrature point
        do ipg = 1, hhoQuad%nbQuadPoints
            coor = hhoCell%barycenter-hhoQuad%points(1:3, ipg)
            qp_coor = hhoQuad%weights(ipg)*coor
            axes(1, 1) = axes(1, 1)+qp_coor(1)*coor(1)
            axes(1, 2) = axes(1, 2)+qp_coor(1)*coor(2)
            axes(1, 3) = axes(1, 3)+qp_coor(1)*coor(3)
            axes(2, 2) = axes(2, 2)+qp_coor(2)*coor(2)
            axes(2, 3) = axes(2, 3)+qp_coor(2)*coor(3)
            axes(3, 3) = axes(3, 3)+qp_coor(3)*coor(3)
        end do
!
! ----- Compute eigenvector
        evalues = 0.d0
        b_n = to_blas_int(hhoCell%ndim)
        b_lda = to_blas_int(3)
        b_lwork = to_blas_int(50)
        call dsyev('V', 'U', b_n, axes, b_lda, &
                   evalues, work, b_lwork, info)
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
        type(HHO_quadrature) :: hhoQuad
        integer(kind=8) :: ipg, idim, i, j, nbpg, is, isn, ie, k
        real(kind=8) :: evalues(3), work(50), axes_3d(3, 3), coor(3)
        real(kind=8) :: xpg(2), poidpg(2), meas_2, weight, coor_pg(3), v1(3)
        real(kind=8) :: bar(3), be, ne(3), xfxe(3), nf(3), f(3), qp_coor(3)
        real(kind=8), parameter :: rac_1div3 = sqrt(1.d0/3.d0)
        blas_int :: b_lda, b_lwork, b_n, info
!
        axes = 0.d0
        axes_3d = 0.d0
!
        if (hhoFace%ndim == 1) then
            coor = hhoFace%coorno(1:3, 2)-hhoFace%coorno(1:3, 1)
            axes(1:3, 1) = coor/norm2(coor)
        else
            if (ASTER_FALSE) then
!
! --- New algo - An efficient method to integrate polynomials over polytopes and curved solids
                nf = hhoFace%normal
!-------- Quadrature of order 3
!            call elraga("SE2", "FPG2", dimp, nbpg, xpg, poidpg)
                nbpg = 2
                poidpg = 1.d0
                xpg(1) = rac_1div3
                xpg(2) = -xpg(1)
!
! ------- Loop on edge
                do ie = 1, hhoFace%nbnodes
                    is = ie
                    isn = ie+1
                    if (isn > hhoFace%nbnodes) then
                        isn = 1
                    end if
                    v1 = (hhoFace%coorno(1:3, isn)-hhoFace%coorno(1:3, is))/2.d0
                    bar = (hhoFace%coorno(1:3, isn)+hhoFace%coorno(1:3, is))/2.d0
                    xfxe = bar-hhoFace%barycenter
                    ! print *, is, isn
                    ! print *, hhoFace%coorno(1:3, isn), hhoFace%coorno(1:3, is)
                    ! print *, nf
                    ! print *, v1
! ----- Compute normal to edge
                    nf = prod_vec(xfxe, v1)
                    ne = prod_vec(nf, v1)
                    if (dot_product(ne, xfxe) < 0.d0) then
                        ne = -ne
                    end if
                    ne = ne/norm2(ne)
! ----- Constant only for plane edge
                    be = dot_product(bar, ne)
                    meas_2 = norm2(v1)
                    ! print *, ne
                    ! print *, be, meas_2
!
! ----- Loop on quadrature point
                    do ipg = 1, nbpg
                        weight = meas_2*poidpg(ipg)
                        coor_pg = bar+xpg(ipg)*v1
                        do i = 1, 3
                            do j = i, 3
                                ! order 0
                                f(1) = hhoFace%barycenter(i)*hhoFace%barycenter(j)
                                ! order 1
                                f(2) = -coor_pg(i)*hhoFace%barycenter(j) &
                                       -coor_pg(j)*hhoFace%barycenter(i)
                                ! order 2
                                f(3) = coor_pg(i)*coor_pg(j)
                                do k = 1, 3
                                    axes_3d(i, j) = axes_3d(i, j)+1.d0/(2.d0+(k-1))*be*weight*f(k)
                                end do
                            end do
                        end do
                    end do
                end do
!
            else
!
! ----- get quadrature
                call hhoQuad%GetQuadFace(hhoFace, 2)
!
! ----- Loop on quadrature point
                do ipg = 1, hhoQuad%nbQuadPoints
                    coor = hhoFace%barycenter-hhoQuad%points(1:3, ipg)
                    ! coor = hhoFace%barycenter
                    qp_coor = hhoQuad%weights(ipg)*coor
                    axes_3d(1, 1) = axes_3d(1, 1)+qp_coor(1)*coor(1)
                    axes_3d(1, 2) = axes_3d(1, 2)+qp_coor(1)*coor(2)
                    axes_3d(1, 3) = axes_3d(1, 3)+qp_coor(1)*coor(3)
                    axes_3d(2, 2) = axes_3d(2, 2)+qp_coor(2)*coor(2)
                    axes_3d(2, 3) = axes_3d(2, 3)+qp_coor(2)*coor(3)
                    axes_3d(3, 3) = axes_3d(3, 3)+qp_coor(3)*coor(3)
                end do
            end if

!
! ----- Compute eigenvector
            evalues = 0.d0
            b_n = to_blas_int(hhoFace%ndim+1)
            b_lda = to_blas_int(3)
            b_lwork = to_blas_int(50)
            call dsyev('V', 'U', b_n, axes_3d, b_lda, &
                       evalues, work, b_lwork, info)
            ASSERT(info == 0)
            ASSERT(minloc(evalues(1:hhoFace%ndim+1), dim=1) == 1)
!
            do idim = 1, hhoFace%ndim
                axes(1:3, idim) = axes_3d(:, idim+1)/norm2(axes_3d(:, idim+1))
            end do
        end if
!
    end function
!
end module
