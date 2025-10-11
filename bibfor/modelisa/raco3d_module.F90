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

module raco3d_module
!
#include "asterf.h"
#include "asterf_types.h"
#include "asterfort/matinv.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/provec.h"

! Global variables (may evolve with time)

    integer(kind=8) :: NB_GAUSS_MAX = 20
    integer(kind=8) :: NB_GREL_MAX = 9
    integer(kind=8) :: NB_NO_CO_MAX = 3
    integer(kind=8) :: NB_NO_3D_MAX = 8
    integer(kind=8) :: NB_NDDL_MAX = 6*4+10*3

! Shell-3d link structures

    type :: PointerContainer
        integer(kind=8), pointer :: iptr(:) => null()
        real(kind=8), pointer :: rptr(:) => null()
    end type PointerContainer
!

contains
!
!

    function segseg_distance(coorseg1, coorseg2)
        ! Compute the shortest distance between two 3D line segments
        real(kind=8), intent(in) :: coorseg1(3, 2), coorseg2(3, 2)
        real(kind=8) :: P1(3), P2(3), Q1(3), Q2(3)
        real(kind=8) :: dist
        real(kind=8) :: u(3), v(3), w(3), a, b, c, d, e
        real(kind=8) :: det, sc, tc
        real(kind=8) :: closestP(3), closestQ(3)
        real(kind=8) :: endpoint_distances(4)
        integer(kind=8) :: i

        do i = 1, 3
            P1(i) = coorseg1(i, 1)
            P2(i) = coorseg1(i, 2)
            Q1(i) = coorseg2(i, 1)
            Q2(i) = coorseg2(i, 2)
        end do

        ! Define vectors along the segments
        u = P2-P1
        v = Q2-Q1
        w = P1-Q1

        ! Dot products
        a = dot_product(u, u)
        b = dot_product(u, v)
        c = dot_product(v, v)
        d = dot_product(u, w)
        e = dot_product(v, w)

        det = a*c-b*b

        ! initialize
        sc = 0.0
        tc = 0.0

        ! Compute the parameters t and s for the closest points
        if (det > 1.0e-8*a) then
            ! Non-parallel case
            sc = (b*e-c*d)/det
            tc = (a*e-b*d)/det
        else if ((-d/a) .ge. 0 .and. (-d/a) .le. 1) then
            sc = -d/a
            tc = 0.0
        else if ((e/c) .ge. 0 .and. (e/c) .le. 1) then
            sc = 0.0
            tc = e/c
        end if

        sc = max(0.0, min(1.0, sc))
        tc = max(0.0, min(1.0, tc))

        ! Compute the closest points
        closestP = P1+sc*u
        closestQ = Q1+tc*v

        ! Compute segment-segment distance
        dist = sqrt(sum((closestP-closestQ)**2))

        ! Compute distances between all endpoints
        endpoint_distances(1) = sqrt(sum((P1-Q1)**2))
        endpoint_distances(2) = sqrt(sum((P1-Q2)**2))
        endpoint_distances(3) = sqrt(sum((P2-Q1)**2))
        endpoint_distances(4) = sqrt(sum((P2-Q2)**2))

        segseg_distance = min(dist, minval(endpoint_distances))

    end function segseg_distance

    function det_jacob(t, s, ddf_co, pg_coor, coor_pt_co, nno_co, epai)
        ! Compute the determinant of the Jacobian at gauss
        ! point of seg(2,3, ..) element
        ! Arguments:
        ! - t, s: Tangent vectors
        ! - ddf_co: second derivatives of shape functions
        ! - pg_coor: Gauss point coordinate
        ! - coor_pt_co: Coordinates of shell points at the boundary
        ! - nno_co: Number of shell nodes
        ! - epai: Thickness of the element
        ! Returns:
        ! - det_jacob: Determinant of the Jacobian

        real(kind=8), intent(in) :: t(3), s(3), coor_pt_co(3, NB_NO_CO_MAX)
        real(kind=8), intent(in) :: ddf_co(:), pg_coor, epai
        integer(kind=8), intent(in) :: nno_co
        !
        real(kind=8) :: magnitude, normal(3)
        real(kind=8) :: vect(3), vect1(3), vect2(3)
        integer(kind=8) :: l, k

        normal = 0.0d0
        vect = 0.0d0
        vect1 = 0.0d0
        vect2 = 0.0d0
        res = 0.0d0

        ! Compute the normal vector as the cross product of s and t
        call provec(s, t, normal)
        magnitude = sqrt(sum(normal(:)**2))
        ! normalize
        normal(:) = normal(:)/magnitude
        !
        ! Compute the vector `vect`
        do l = 1, 3
            do k = 1, nno_co
                vect(l) = vect(l)+ddf_co(k)*coor_pt_co(l, k)
            end do
        end do

        ! Compute vect1 as the cross product of s and vect
        call provec(s, vect, vect1)

        ! Project vect1 orthogonally to the plane defined by the normal vector
        vect1(:) = (vect1(:)-dot_product(vect1(:), normal(:))*normal(:))/magnitude

        ! Compute vect2 as a combination of t and a thickness term involving vect1
        vect2(:) = t(:)+0.5*epai*pg_coor*vect1(:)

        ! Compute the determinant of the Jacobian
        det_jacob = 0.5*epai*sqrt(sum(vect2(:)**2))

    end function det_jacob

    function find_parametric_coordinates(p, typma3d, coor_pt_3d) result(res)
        implicit none
        ! Inputs
        character(len=8), intent(in) :: typma3d
        real(kind=8), intent(in) :: p(3)
        real(kind=8), intent(in) :: coor_pt_3d(:, :)

        ! Outputs
        real(kind=8) :: res(3)
        ! Others
        real(kind=8) :: xi(2), dxi(2), ff(8), df(3, 8)
        integer(kind=8) :: max_iter, iter
        real(kind=8) :: F(3), jac(3, 2), det, residual(2)
        real(kind=8) ::  tol
        real(kind=8), dimension(2, 2) :: JtJ
        real(kind=8), dimension(2, 3) :: Jt
        real(kind=8), dimension(2, 2) :: JtJ_inv
        aster_logical :: converged

        ! Initialize
        tol = 1.0e-6
        max_iter = 30
        res = 0.d0
        res(3) = -1.d0
        ff = 0.d0
        df = 0.d0
        xi = 0.d0
        dxi = 0.d0
        F = 0.d0
        converged = .false.
        iter = 0

        do while (.not. converged .and. iter < max_iter)

            call elrfvf(typma3d, xi, ff)
            call elrfdf(typma3d, xi, df)

            F(1) = sum(ff*coor_pt_3d(1, :))-p(1)
            F(2) = sum(ff*coor_pt_3d(2, :))-p(2)
            F(3) = sum(ff*coor_pt_3d(3, :))-p(3)

            ! Jocobian
            jac(1, 1) = sum(df(1, :)*coor_pt_3d(1, :))
            jac(1, 2) = sum(df(2, :)*coor_pt_3d(1, :))
            jac(2, 1) = sum(df(1, :)*coor_pt_3d(2, :))
            jac(2, 2) = sum(df(2, :)*coor_pt_3d(2, :))
            jac(3, 1) = sum(df(1, :)*coor_pt_3d(3, :))
            jac(3, 2) = sum(df(2, :)*coor_pt_3d(3, :))

            ! Compute J^T (transpose of J)
            Jt = transpose(jac)

            ! Compute J^T * J (2x2 matrix)
            JtJ = matmul(Jt, jac)

            ! Compute the inverse of J^T * J (2x2 matrix)
            call matinv('S', 2, JtJ, JtJ_inv, det)

            ! residual
            residual = -matmul(Jt, F)

            ! solve
            dxi = matmul(JtJ_inv, residual)

            ! update
            xi = xi+dxi

            if (maxval(abs(dxi)) < tol) then
                converged = .true.
            end if

            iter = iter+1

        end do

        ! Last check before the return
        if (converged) then
            if (((typma3d(1:3) .eq. "TR3" .or. typma3d(1:3) .eq. "TR6") &
                 .and. xi(1) >= 0.0d0 .and. xi(1) <= 1.0d0 .and. &
                 xi(2) >= 0.0d0 .and. xi(2) <= 1.0d0) .or. &
                ((typma3d(1:3) .eq. "QU4" .or. typma3d(1:3) .eq. "QU8") &
                 .and. xi(1) >= -1.0d0 .and. xi(1) <= 1.0d0 .and. &
                 xi(2) >= -1.0d0 .and. xi(2) <= 1.0d0)) then
                ! Return the coordinates
                res(1) = xi(1)
                res(2) = xi(2)
                res(3) = 1.d0
            end if
        end if

    end function find_parametric_coordinates

!
end module raco3d_module
