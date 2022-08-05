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
use contact_type
use contact_algebra_module
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
#include "asterfort/mm2onf.h"
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
!
    public :: projQpSl2Ma, shapeFuncDisp, shapeFuncLagr, evalPoly
    public :: diameter, testLagrC, gapEval, testLagrF
    public :: speedEval, thresEval
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine projQpSl2Ma(geom, coor_qp_sl, proj_tole, &
                            coor_qp_ma, gap, &
                            tau_slav, norm_slav, tau_mast, norm_mast)
!
    implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_sl(2), proj_tole
        real(kind=8), intent(out) :: coor_qp_ma(2), gap
        real(kind=8), intent(out) :: tau_slav(3,2), tau_mast(3,2)
        real(kind=8), intent(out) :: norm_slav(3), norm_mast(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Project slave point to master space and compute quantity
!
! --------------------------------------------------------------------------------------------------
!
        integer :: iret, elem_mast_line_nbnode
        real(kind=8) :: coor_qp_sl_re(3), tau1_mast(3), tau2_mast(3), coor_qp_ma_re(3)
        real(kind=8) :: ksi1_line, ksi2_line
        character(len=8) :: elem_mast_line_code
!
        coor_qp_ma = 0.d0
        norm_slav = 0.d0
        norm_mast = 0.d0
        tau_slav = 0.d0
        tau_mast = 0.d0
        gap = 0.d0
!
! ------ Compute outward slave normal (pairing configuration)
!
        call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_pair, &
                    coor_qp_sl(1), coor_qp_sl(2), norm_slav, tau_slav(1:3,1), tau_slav(1:3,2))
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
!
! ----- Initialize with linearized cell
!
        if(iret == 1) then
            if(geom%elem_mast_code(1:2) == "SE") then
                elem_mast_line_code = "SE2"
                elem_mast_line_nbnode = 2
            elseif(geom%elem_mast_code(1:2) == "TR") then
                elem_mast_line_code = "TR3"
                elem_mast_line_nbnode = 3
            elseif(geom%elem_mast_code(1:2) == "QU") then
                elem_mast_line_code = "QU4"
                elem_mast_line_nbnode = 4
            else
                ASSERT(ASTER_FALSE)
            end if

            call mmnewd(elem_mast_line_code, elem_mast_line_nbnode, geom%elem_dime, &
                        geom%coor_mast_curr, coor_qp_sl_re, 75, proj_tole, norm_slav, &
                        ksi1_line, ksi2_line, tau1_mast, tau2_mast, iret)
            ASSERT(iret==0)
!
            call mmnewd(geom%elem_mast_code, geom%nb_node_mast, geom%elem_dime, &
                        geom%coor_mast_curr, coor_qp_sl_re, 75, proj_tole, norm_slav, &
                        coor_qp_ma(1), coor_qp_ma(2), tau1_mast, tau2_mast, iret, &
                        ksi1_init=ksi1_line, ksi2_init= ksi2_line)
            ASSERT(iret==0)
!
        end if
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
                    coor_qp_ma(1), coor_qp_ma(2), norm_mast, tau_mast(1:3,1), tau_mast(1:3,2))
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
    subroutine thresEval(param, l_cont_qp, projRmVal, thres_qp, l_fric_qp)
!
    implicit none
!
        type(ContactParameters), intent(in) :: param
        aster_logical, intent(in) :: l_cont_qp
        real(kind=8), intent(in) :: projRmVal
        aster_logical, intent(out) :: l_fric_qp
        real(kind=8), intent(out) :: thres_qp
!
! --------------------------------------------------------------------------------------------------
!
!   Compute threshold for friction
!
! --------------------------------------------------------------------------------------------------
!
        thres_qp = 0.d0
        l_fric_qp = ASTER_FALSE
!
! ----- Define threshold for friction
        if(param%l_fric) then
            if(param%type_fric == FRIC_TYPE_TRES) then
                l_fric_qp = ASTER_TRUE
                thres_qp = param%threshold_given
            elseif(param%type_fric == FRIC_TYPE_NONE) then
                l_fric_qp = ASTER_FALSE
                thres_qp = 0.d0
            else
                l_fric_qp = l_cont_qp
                if(l_cont_qp) then
                    if (param%type_fric == FRIC_TYPE_COUL) then
                        thres_qp = - param%threshold_given * projRmVal
                    else
                        thres_qp = THRES_STICK
                    end if
                else
                    thres_qp = 0.d0
                end if
            end if
        end if
!
    end subroutine
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
    subroutine shapeFuncDisp(elem_dime, elem_nbnode, elem_code, coor_qp, &
                             shape_, dshape_, ddshape_)
!
    implicit none
!
        integer, intent(in) :: elem_dime
        integer, intent(in) :: elem_nbnode
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: coor_qp(2)
        real(kind=8), intent(out), optional :: shape_(9)
        real(kind=8), intent(out), optional  :: dshape_(2,9)
        real(kind=8), intent(out), optional  :: ddshape_(3,9)
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
!
        if(present(ddshape_)) then
            call mm2onf(elem_dime, elem_nbnode, elem_code, coor_qp(1), coor_qp(2), ddshape_)
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
end module
