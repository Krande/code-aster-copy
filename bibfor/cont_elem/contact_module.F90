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
module contact_module
!
!
implicit none
!
private
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmnewd.h"
#include "asterfort/mmnonf.h"
#include "asterfort/mmnorm.h"
#include "asterfort/reerel.h"
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
! Contact - generic
!
! Generic method for contact
!
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
!
type ContactParameters
    !! Contact parameters
    integer                             :: algo_cont = 0
    integer                             :: type_cont = 0
    real(kind=8)                        :: vari_cont = 0.d0
    real(kind=8), dimension(4)          :: coef_cont = 0.d0

    !! Friction paramaters
    aster_logical                       :: l_fric = ASTER_FALSE
    integer                             :: algo_fric = 0
    integer                             :: type_fric = 0
    real(kind=8), dimension(4)          :: coef_fric = 0.d0
    real(kind=8), dimension(4)          :: threshold = 0.d0
    real(kind=8)                        :: threshold_given = 0.d0

    !! Other
    real(kind=8)                        :: proj_tole = 0.d0
end type
!
type ContactGeom
    !! Slave side parameters
    integer                             :: nb_node_slav = 0
    character(len=8)                    :: elem_slav_code = " "
    real(kind=8), dimension(3,9)        :: coor_slav_init = 0.d0
    real(kind=8), dimension(3,9)        :: coor_slav_prev = 0.d0
    real(kind=8), dimension(3,9)        :: coor_slav_curr = 0.d0
    real(kind=8), dimension(3,9)        :: coor_slav_pair = 0.d0
    real(kind=8), dimension(3,9)        :: depl_slav_curr = 0.d0
    real(kind=8), dimension(4)          :: lagc_slav_curr = 0.d0
    real(kind=8), dimension(2,4)        :: lagf_slav_curr = 0.d0
    integer                             :: nb_lagr_c      = 0
    integer, dimension(9)               :: indi_lagc      = 0

    !! Master side paramaters
    integer                             :: nb_node_mast = 0
    character(len=8)                    :: elem_mast_code = " "
    real(kind=8), dimension(3,9)        :: coor_mast_init = 0.d0
    real(kind=8), dimension(3,9)        :: coor_mast_prev = 0.d0
    real(kind=8), dimension(3,9)        :: coor_mast_curr = 0.d0
    real(kind=8), dimension(3,9)        :: coor_mast_pair = 0.d0
    real(kind=8), dimension(3,9)        :: depl_mast_curr = 0.d0

    !! Time
    real(kind=8) :: time_prev = 0.d0, time_curr = 0.d0

    !! Other
    integer                             :: elem_dime = 0
    aster_logical                       :: l_axis = ASTER_FALSE
    integer                             :: nb_dofs = 0
end type
!
!===================================================================================================
!
!===================================================================================================
!
    public :: ContactParameters, ContactGeom
    public :: projQpSl2Ma, Heaviside, shapeFuncDisp, shapeFuncLagr, evalPoly
    public :: dGap_du, d2Gap_du2, diameter, testLagrC, gapEval, testLagrF
    public :: projBs, projRm, projTn, jump_tang, jump_norm, dNs_du, otimes, Iden3
    public :: dProjBs_dx, dprojBs_ds, dprojBs_dn, speedEval
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine projQpSl2Ma(geom, coor_qp_sl, proj_tole, &
                            coor_qp_ma, gap, &
                            tau_1_slav, tau_2_slav, norm_slav, norm_mast)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_sl(2), proj_tole
        real(kind=8), intent(out) :: coor_qp_ma(2), gap
        real(kind=8), intent(out) :: tau_1_slav(3), tau_2_slav(3)
        real(kind=8), intent(out) :: norm_slav(3), norm_mast(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Project slave point to master space and compute quantity
!
! --------------------------------------------------------------------------------------------------
!
        integer :: iret
        real(kind=8) :: coor_qp_sl_re(3), tau1_mast(3), tau2_mast(3), coor_qp_ma_re(3)
!
        coor_qp_ma = 0.d0
        norm_slav = 0.d0
        norm_mast = 0.d0
        tau_1_slav = 0.d0
        tau_2_slav = 0.d0
        gap = 0.d0
!
! ------ Compute outward slave normal (pairing configuration)
!
        call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_pair, &
                    coor_qp_sl(1), coor_qp_sl(2), norm_slav)
!
! ----- Return in real slave space (current configuration)
!
        coor_qp_sl_re = 0.d0
        call reerel(geom%elem_slav_code, geom%nb_node_slav, 3, geom%coor_slav_curr, coor_qp_sl, &
                    coor_qp_sl_re)
!
! ----- Projection of node on master cell (master parametric space)
!
        call mmnewd(geom%elem_mast_code, geom%nb_node_mast, geom%elem_dime, geom%coor_mast_curr,&
                    coor_qp_sl_re, 75, proj_tole, norm_slav, &
                    coor_qp_ma(1), coor_qp_ma(2), tau1_mast, tau2_mast, &
                    iret)
        ASSERT(iret == 0)
!
! ----- Return in real master space
!
        coor_qp_ma_re = 0.d0
        call reerel(geom%elem_mast_code, geom%nb_node_mast, 3, geom%coor_mast_curr, coor_qp_ma, &
                    coor_qp_ma_re)
!
! ------ Compute outward master normal (pairing configuration)
!

        call apnorm(geom%nb_node_mast, geom%elem_mast_code, geom%elem_dime, geom%coor_mast_pair, &
                    coor_qp_ma(1), coor_qp_ma(2), norm_mast)
!
! ----- Compute gap for raytracing gap = -(x^s - x^m).n^s
!
        gap = gapEval(coor_qp_sl_re, coor_qp_ma_re, norm_slav)

        ! print*, "COOR_SL: ", geom%coor_slav_curr(1,1:2)
        ! print*, "COOR_MA: ", geom%coor_mast_curr(1,1:2)
        ! print*, "NORM_SL: ", norm_slav
        ! print*, "NORM_MA: ", norm_mast
        ! print*, "COOR_QP: ", coor_qp_sl
        ! print*, "COOR_QP_RE: ", coor_qp_sl_re
        ! print*, "COOR_PJ: ", coor_qp_ma
        ! print*, "COOR_PJ_RE: ", coor_qp_ma_re
        ! print*, "GAP: ", gap
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function gapEval(slav_pt, mast_pt, norm_slav)
!
    implicit none
!
        real(kind=8), intent(in) :: slav_pt(3), mast_pt(3), norm_slav(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute gap for raytracing gap = -(x^s - x^m).n^s
!
! --------------------------------------------------------------------------------------------------
!
        gapEval = -dot_product(slav_pt-mast_pt, norm_slav)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function speedEval(geom, coor_qp_slav, coor_qp_mast, gap)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_slav(2), coor_qp_mast(2), gap
        real(kind=8) :: speedEval(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute speed for raytracing v = -(x^s(t_prev) - x^m(t_prev) + gap * n^s(t_prev))
!                                     /(t_curr - t_prev)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: coor_qp_sl_prev(3), coor_qp_ma_prev(3), norm_slav_prev(3)
!
! ----- Return in real slave space
!
        coor_qp_sl_prev = 0.d0
        call reerel(geom%elem_slav_code, geom%nb_node_slav, 3, geom%coor_slav_prev, coor_qp_slav, &
                    coor_qp_sl_prev)
!
! ----- Return in real master space
!
        coor_qp_ma_prev = 0.d0
        call reerel(geom%elem_mast_code, geom%nb_node_mast, 3, geom%coor_mast_prev, coor_qp_mast, &
                    coor_qp_ma_prev)
!
! ------ Compute outward slave normal
!
        call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_prev, &
                    coor_qp_slav(1), coor_qp_slav(2), norm_slav_prev)
!
        speedEval = -(coor_qp_sl_prev - coor_qp_ma_prev + gap * norm_slav_prev) &
                    / (geom%time_curr - geom%time_prev)
!
    end function
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
        Iden3 = 0
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
    subroutine shapeFuncDisp(elem_dime, elem_nbnode, elem_code, coor_qp, &
                             shape_, dshape_)
!
    implicit none
!
        integer, intent(in) :: elem_dime
        integer, intent(in) :: elem_nbnode
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: coor_qp(2)
        real(kind=8), intent(out), optional :: shape_(9)
        real(kind=8), intent(out), optional  :: dshape_(2,9)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate shape function and derivative for displacement
!
! --------------------------------------------------------------------------------------------------
!
        if(present(shape_)) then
            call mmnonf(elem_dime, elem_nbnode, elem_code, coor_qp(1), coor_qp(2), shape_)
        end if
!
        if(present(dshape_)) then
            call mmdonf(elem_dime, elem_nbnode, elem_code, coor_qp(1), coor_qp(2), dshape_)
        end if
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine shapeFuncLagr(elem_dime, elem_code, coor_qp, &
                             shape_)
!
    implicit none
!
        integer, intent(in) :: elem_dime
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: coor_qp(2)
        real(kind=8), intent(out) :: shape_(4)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate shape function and derivative (P1 for Lagrange)
!
! --------------------------------------------------------------------------------------------------
!
        character(len=8) :: elem_code_lagr
        integer :: elem_nbnode_lagr
        real(kind=8) :: ff(9)
!
        if(elem_code == "SE2") then
            elem_code_lagr = "SE2"
            elem_nbnode_lagr = 2
        elseif(elem_code == "SE3") then
                elem_code_lagr = "SE3"
                elem_nbnode_lagr = 3
        elseif(elem_code(1:2) == "TR") then
            elem_code_lagr = "TR3"
            elem_nbnode_lagr = 3
        elseif(elem_code(1:2) == "QU") then
            elem_code_lagr = "QU4"
            elem_nbnode_lagr = 4
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call mmnonf(elem_dime, elem_nbnode_lagr, elem_code_lagr, &
                        coor_qp(1), coor_qp(2), ff)
        shape_(1:4) = ff(1:4)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function evalPoly(nb_node, shape, coeff_node)
!
    implicit none
!
        integer, intent(in) :: nb_node
        real(kind=8), intent(in) :: coeff_node(*)
        real(kind=8), intent(in) :: shape(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate polynome via linear combinaison
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node
!
        evalPoly = 0.d0
        do i_node= 1, nb_node
            evalPoly = evalPoly + coeff_node(i_node) * shape(i_node)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function diameter(nb_node, nodes_coor)
!
    implicit none
!
        integer, intent(in) :: nb_node
        real(kind=8), intent(in) :: nodes_coor(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate diameter of cell
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node_i, i_node_j
        real(kind=8) :: length
!
        diameter = 0.d0
        do i_node_i= 1, nb_node
            do i_node_j= i_node_i + 1, nb_node
                length = norm2(nodes_coor(1:3, i_node_i) - nodes_coor(1:3, i_node_j))
                diameter = max(diameter, length)
            end do
        end do
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
!===================================================================================================
!
!===================================================================================================
!
    subroutine testLagrC(geom, func_lagr, mu_c)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_lagr(4)
        real(kind=8), intent(out) :: mu_c(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate test Lagrangian function
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node, index
!
        mu_c = 0.d0
        index = 0
!
! --- Slave side
!
        do i_node = 1, geom%nb_lagr_c
            ASSERT(geom%indi_lagc(i_node) > 0)
            index = index + geom%elem_dime
            mu_c(index+1) = func_lagr(i_node)
            index = index + geom%indi_lagc(i_node)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine testLagrF(geom, func_lagr, mu_f)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_lagr(4)
        real(kind=8), intent(out) :: mu_f(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate test Lagrangian function
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_node, index
!
        mu_f = 0.d0
        index = 0
!
! --- Slave side
!
        do i_node = 1, geom%nb_lagr_c
            ASSERT(geom%indi_lagc(i_node) > 0)
            index = index + geom%elem_dime
            mu_f(index+2, 1) = func_lagr(i_node)
            if(geom%elem_dime == 3) then
                mu_f(index+3, 2) = func_lagr(i_node)
            end if
            index = index + geom%indi_lagc(i_node)
        end do
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
