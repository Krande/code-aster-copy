! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
module contact_algebra_module
!
use contact_type
!
implicit none
!
private
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/subac1.h"
#include "asterfort/subaco.h"
#include "asterfort/subacv.h"
#include "asterfort/sumetr.h"
#include "blas/dgemv.h"
#include "blas/dger.h"
#include "contact_module.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Linear algebra and basic derivatives
!
! --------------------------------------------------------------------------------------------------
!
    public :: Heaviside, jump_norm, jump_tang
    public :: dGap_du, d2Gap_du2
    public :: projBs, projRm, projTn, dNs_du, otimes, Iden3
    public :: dProjBs_dx, dprojBs_ds, dprojBs_dn
    public :: metricTensor, invMetricTensor, cmpTang
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function projRm(x)
!
    implicit none
!
        real(kind=8), intent(in) :: x
!
! --------------------------------------------------------------------------------------------------
!
!   Project x on R^- = min(x, 0)
!
! --------------------------------------------------------------------------------------------------
!
        if(x <= TOLE_BORNE) then
            projRm = x
        else
            projRm = 0.d0
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function Heaviside(x)
!
    implicit none
!
        real(kind=8), intent(in) :: x
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate Heaviside function
!
! --------------------------------------------------------------------------------------------------
!
        if(x >= -(TOLE_BORNE)) then
            Heaviside = 1.d0
        else
            Heaviside = 0.d0
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function otimes(v1, v2)
!
    implicit none
!
        real(kind=8), intent(in) :: v1(3), v2(3)
        real(kind=8) :: otimes(3,3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate otimes product: (v1 \otimes v2 )_{i,j} = v1_{i} * v2_{j}
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i, j
!
        do j = 1, 3
            do i= 1, 3
                otimes(i,j) = v1(i) * v2(j)
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function Iden3()
!
    implicit none
!
        real(kind=8) :: Iden3(3,3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate identity matrix
!
! --------------------------------------------------------------------------------------------------
!
        Iden3 = 0.d0
        Iden3(1,1) = 1.d0
        Iden3(2,2) = 1.d0
        Iden3(3,3) = 1.d0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function projTn(normal)
!
    implicit none
!
        real(kind=8), intent(in) :: normal(3)
        real(kind=8) :: projTn(3,3)
!
! --------------------------------------------------------------------------------------------------
!
!   Projection operator onto the tangent plane corresponding to normal vector n
!   Tn = Id - n \otimes n
!
! --------------------------------------------------------------------------------------------------
!
        projTn = Iden3() - otimes(normal, normal)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function metricTensor(tau)
!
    implicit none
!
        real(kind=8), intent(in) :: tau(3,2)
        real(kind=8) :: metricTensor(2,2)
!
! --------------------------------------------------------------------------------------------------
!
!   Metric tensor: m[i,j] = tau_i . tau_j
!
! --------------------------------------------------------------------------------------------------
!
        metricTensor = matmul(transpose(tau), tau)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function invMetricTensor(geom, metricTens)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: metricTens(2,2)
        real(kind=8) :: invMetricTensor(2,2)
!
! --------------------------------------------------------------------------------------------------
!
!   Metric tensor: m[i,j] = tau_i . tau_j
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: det
!
        invMetricTensor = 0.d0
!
        if(geom%elem_dime == 3) then
            det = metricTens(1,1) * metricTens(2,2) - metricTens(1,2) * metricTens(2,1)
            invMetricTensor(1,1) = metricTens(2,2) / det
            invMetricTensor(1,2) = -metricTens(2,1) / det
            invMetricTensor(2,1) = invMetricTensor(1,2)
            invMetricTensor(2,2) = metricTens(1,1) / det
        else
            invMetricTensor(1,1) = 1.d0 / metricTens(1,1)
         end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function cmpTang(invMetricTens, tau, vec)
!
    implicit none
!
        real(kind=8), intent(in) :: invMetricTens(2,2)
        real(kind=8), intent(in) :: tau(3,2), vec(3)
        real(kind=8) :: cmpTang(2)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute coefficient in tangential basis (tau_1, tau_2)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: rhs(2)
!
        rhs = matmul(transpose(tau), vec)
!
        cmpTang = matmul(invMetricTens, rhs)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function projBs(param, x, s, norm_slav)
!
    implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), norm_slav(3)
        real(kind=8), intent(in) :: s
        real(kind=8) :: projBs(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Project x on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3,3), xT(3)
!
        projBs = 0
        if(param%type_fric .ne. FRIC_TYPE_NONE) then
            if(abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
!
                if(norm_xT <= s) then
                    projBs = xT
                else
                    projBs = s * xT / norm_xT
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojBs_dx(param, x, s, norm_slav)
!
    implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3)
        real(kind=8), intent(in) :: s, norm_slav(3)
        real(kind=8) :: dprojBs_dx(3,3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative along x of the projection on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3,3), xT_uni(3), xT(3)
!
        dprojBs_dx = 0
        if(param%type_fric .ne. FRIC_TYPE_NONE) then
            if(abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
!
                if(norm_xT <= s) then
                    dprojBs_dx = Tn
                else
                    xT_uni = xT / norm_xT
                    dprojBs_dx = s * norm_xT * (Tn - otimes(xT_uni, xT_uni))
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojBs_ds(param, x, s, norm_slav)
!
    implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3), norm_slav(3)
        real(kind=8), intent(in) :: s
        real(kind=8) :: dprojBs_ds(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative along s of the projection on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3,3), xT(3)
!
        dprojBs_ds = 0
        if(param%type_fric .ne. FRIC_TYPE_NONE) then
            if(abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
!
                if(norm_xT > s) then
                    dprojBs_ds = xT / norm_xT
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function dprojBs_dn(param, x, s, norm_slav)
!
    implicit none
!
        type(ContactParameters), intent(in) :: param
        real(kind=8), intent(in) :: x(3)
        real(kind=8), intent(in) :: s, norm_slav(3)
        real(kind=8) :: dprojBs_dn(3,3)
!
! --------------------------------------------------------------------------------------------------
!
!   Derivative along the normal of the projection on the sphere center to zero of radius s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm_xT, Tn(3,3), xT_uni(3), x_n, xT(3)
!
        dprojBs_dn = 0
        if(param%type_fric .ne. FRIC_TYPE_NONE) then
            if(abs(s) > r8prem()) then
                Tn = projTn(norm_slav)
                xT = matmul(Tn, x)
                norm_xT = norm2(xT)
                x_n = dot_product(x, norm_slav)
!
                if(norm_xT <= s) then
                    dprojBs_dn = -x_n * Tn - otimes(norm_slav, xT)
                else
                    xT_uni = xT / norm_xT
                    dprojBs_dn = -s/norm_xT*(x_n*(Tn-otimes(xT_uni,xT_uni)) + otimes(norm_slav, xT))
                end if
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dNs_du(geom, func_slav, dfunc_slav, dNs)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_slav(2,9), func_slav(9)
        real(kind=8), intent(out) :: dNs(MAX_LAGA_DOFS,3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of slave normal using differential geometry
!   D n^s(u^s)[v^s] = sum_{i=1, d-1}[dv_di.g^i * n^s - dv_di.n^s * g^i]
!   with g^i - controvariant basis
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: cova(3, 3), metr(2, 2), jac, cnva(3, 2), acv(2, 2)
        real(kind=8) :: t1(3), t2(3), t3(3), norm_slav(3), ta(3)
        integer :: i_node_slav, i_elem_dime, index
!
! ----- Covariant basis
!
        if(geom%elem_dime == 3) then
            call subaco(geom%nb_node_slav, dfunc_slav, geom%coor_slav_pair, cova)
        else
            call subac1(geom%l_axis, geom%nb_node_slav, func_slav, dfunc_slav(1,:), &
                        geom%coor_slav_pair(1:2,:), cova)
        end if
!
! ----- Metric tensor
!
        call sumetr(cova, metr, jac)
!
! ----- Contra-variant basis
!
        call subacv(cova, metr, jac, cnva, acv)
!
        norm_slav(1:3) = cova(1:3,3)
!
! ----- Tangent matrix
!
        index = 0
        dNs = 0.d0
        ta = 0.d0
        do i_node_slav = 1, geom%nb_node_slav
            do i_elem_dime = 1, geom%elem_dime
                index = index + 1
                t1 = (dfunc_slav(1,i_node_slav)*cnva(i_elem_dime,1) + &
                      dfunc_slav(2,i_node_slav)*cnva(i_elem_dime,2)) * norm_slav
                t2 = dfunc_slav(1,i_node_slav)*norm_slav(i_elem_dime) * cnva(1:3,1)
                t3 = dfunc_slav(2,i_node_slav)*norm_slav(i_elem_dime) * cnva(1:3,2)
                if(geom%l_axis .and. i_elem_dime == 1) then
                    ta = func_slav(i_node_slav)*cnva(3,2)*norm_slav
                end if
                dNs(index, 1:3) = dNs(index, 1:3) + (t1 - t2 - t3 + ta)

            enddo
            index = index + geom%indi_lagc(i_node_slav)
        enddo
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dGap_du(geom, func_slav, dfunc_slav, func_mast, &
                       norm_slav, norm_mast, gap, dGap)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_slav(9), func_mast(9), dfunc_slav(2,9), gap
        real(kind=8), intent(in) :: norm_slav(3), norm_mast(3)
        real(kind=8), intent(out) :: dGap(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate first derivative of gap for raytracing
!   D gap(u)[v] = -(v^s - v^m + gap*Dn^s[v^s]).n^m/(n^m.n^s)
!   ou
!   D gap(u)[v] = -(v^s - v^m + dy_de * De[v]).n^s
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: norm(3), dNs(MAX_LAGA_DOFS,3)
!
        norm = norm_mast / dot_product(norm_mast, norm_slav)
!
! --- normal jump : -(v^s - v^m ).norm
!
        call jump_norm(geom, func_slav, func_mast, norm, dGap)
!
! --- Compute term: -g*Dn^s.norm (= 0 if norm = n^s)
!
        call dNs_du(geom, func_slav, dfunc_slav, dNs)
        call dgemv('N', geom%nb_dofs, geom%elem_dime, -gap, dNs, MAX_LAGA_DOFS, &
                        norm, 1, 1.d0, dGap, 1)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine d2Gap_du2(geom, func_slav, dfunc_slav, func_mast, &
                         norm_slav, norm_mast, gap, d2Gap)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_slav(9), func_mast(9), dfunc_slav(2,9), gap
        real(kind=8), intent(in) :: norm_slav(3), norm_mast(3)
        real(kind=8), intent(out) :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate second derivative of gap
!   D^2 gap(u)[v,w] =
!   - ( D gap(u)[v]* D n^s[w^s] + D gap(u)[w]* D n^s[v^s] + gap(u) D^2 n^s[v^s, w^s] ).n^m/(n^m.n^s)
!   + ( D dy/de[v] * De[w] + D dy/de[w] * De[v] + d^2y/de^2 * De[v]*De[w] ).n^m/(n^m.n^s)
!
!   ou
!   D^2 gap(u)[v,w] =
!   - gap(u) D^2 n^s[v^s, w^s].n^s
!   + ( D dy/de[v] * De[w] + D dy/de[w] * De[v] + d^2y/de^2 * De[v]*De[w] + dy_de * D^2 e[v,w]).n^s
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: vers1 = ASTER_TRUE
        real(kind=8) :: norm(3), dNs(MAX_LAGA_DOFS,3), dGap(MAX_LAGA_DOFS), dNs_n(MAX_LAGA_DOFS)
!
        d2Gap = 0.d0
!
! --- Term: n^m/(n^m.n^s)
!
        if(vers1) then
            norm = norm_mast / dot_product(norm_mast, norm_slav)
        else
            norm = norm_slav
        end if
!
        if(vers1) then
!
! --- Term: D gap(u)[v]
!
            call dGap_du(geom, func_slav, dfunc_slav, func_mast, norm_slav, norm_mast, gap, dGap)
!
! --- Term: D( n^s[v^s]).n^m/(n^m.n^s)
!
            dNs_n = 0.d0
            call dNs_du(geom, func_slav, dfunc_slav, dNs)
            call dgemv('N', geom%nb_dofs, geom%elem_dime, 1.d0, dNs, MAX_LAGA_DOFS, &
                    norm, 1, 1.d0, dNs_n, 1)
!
! --- Term: -(D gap(u)[v]* D n^s[w^s] + D gap(u)[w]* D n^s[v^s]).n^m/(n^m.n^s)
!
            call dger(geom%nb_dofs, geom%nb_dofs, -1.d0, dGap, 1, dNs_n, 1, d2Gap, MAX_LAGA_DOFS)
            call dger(geom%nb_dofs, geom%nb_dofs, -1.d0, dNs_n, 1, dGap, 1, d2Gap, MAX_LAGA_DOFS)
        else
!
! --- Term: dy_de * D^2 e[v,w]
!
        end if
!
! --- Term: D^2 n^s[v^s, w^s]
!

!
! --- Term: D dy/de[v]
!

!
! --- Term: De[v]
!

!
! --- Il manque les autres termes mais pas facile à calculer.
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine jump_tang(geom, func_slav, func_mast, &
                       norm_slav, jump_t)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_slav(9), func_mast(9)
        real(kind=8), intent(in) :: norm_slav(3)
        real(kind=8), intent(out) :: jump_t(MAX_LAGA_DOFS, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate tangential jum of test function:
!   jump_tang = (v^s - v^m).tau_i_slav
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node, i_dim, index, j_dim
        real(kind=8) :: Tn(3,3)
!
        jump_t = 0.d0
        index = 0
        Tn = projTn(norm_slav)
!
! --- Slave side
!
        do i_node = 1, geom%nb_node_slav
            do i_dim = 1, geom%elem_dime
                index = index + 1
                do j_dim = 1, geom%elem_dime
                    jump_t(index, j_dim) = func_slav(i_node) * Tn(j_dim, i_dim)
                end do
            end do
!
            index = index + geom%indi_lagc(i_node)
        end do
!
! --- Master side
!
        do i_node = 1, geom%nb_node_mast
            do i_dim = 1, geom%elem_dime
                index = index + 1
                do j_dim = 1, geom%elem_dime
                    jump_t(index, j_dim) = -func_mast(i_node) * Tn(j_dim, i_dim)
                end do
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine jump_norm(geom, func_slav, func_mast, norm_slav, jump_n)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_slav(9), func_mast(9)
        real(kind=8), intent(in) :: norm_slav(3)
        real(kind=8), intent(out) :: jump_n(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate normal jump of test function:
!   jump_norm = -(v^s - v^m).norm^s
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node, i_dim, index
!
        jump_n = 0.d0
        index = 0
!
! --- Slave side
!
        do i_node = 1, geom%nb_node_slav
            do i_dim = 1, geom%elem_dime
                index = index + 1
                jump_n(index) = -func_slav(i_node) * norm_slav(i_dim)
            end do
!
            index = index + geom%indi_lagc(i_node)
        end do
!
! --- Master side
!
        do i_node = 1, geom%nb_node_mast
            do i_dim = 1, geom%elem_dime
                index = index + 1
                jump_n(index) = func_mast(i_node) * norm_slav(i_dim)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
