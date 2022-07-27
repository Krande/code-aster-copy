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
subroutine laElemCont(parameters, geom, coor_qp_sl, hF, &
                    lagr_c, gap, gamma_c, projRmVal, l_cont_qp,&
                    lagr_f, vT, gamma_f, projBsVal, l_fric_qp, &
                    dGap, d2Gap, mu_c, mu_f, jump_t)
!
use contact_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "contact_module.h"
!
type(ContactParameters), intent(in) :: parameters
type(ContactGeom), intent(in) :: geom
real(kind=8), intent(in) :: coor_qp_sl(2), hF
real(kind=8), intent(out) :: lagr_c, gap, gamma_c, projRmVal
real(kind=8), intent(out) :: lagr_f(2), vT(2), gamma_f, projBsVal(2)
aster_logical, intent(out) :: l_cont_qp, l_fric_qp
real(kind=8), intent(out), optional :: dGap(MAX_LAGA_DOFS)
real(kind=8), intent(out), optional :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
real(kind=8), intent(out), optional :: mu_c(MAX_LAGA_DOFS)
real(kind=8), intent(out), optional :: mu_f(MAX_LAGA_DOFS,2)
real(kind=8), intent(out), optional :: jump_t(MAX_LAGA_DOFS,2)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute elementary contact quantities
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_dime        : dimension of elements
! In  l_axis           : .true. for axisymmetric element
! In  nb_lagr          : total number of Lagrangian dof on contact element
! In  indi_lagc        : node where Lagrangian dof is present (1) or not (0)
! In  lagr_c           : value of contact lagrangian
! In  l_norm_smooth    : indicator for normals smoothing
! In  nb_node_slav     : number of nodes of for slave side from contact element
! In  elem_slav_code   : code element for slave side from contact element
! In  nb_node_mast     : number of nodes of for master side from contact element
! In  elem_mast_code   : code element for master side from contact element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: shape_func_sl(9), shape_func_ma(9), shape_func_lagr(4), dshape_func_sl(2,9)
    real(kind=8) :: norm_slav(3), norm_mast(3), H, coor_qp_ma(2), lagrc_gap, lagrf_vT(2)
    real(kind=8) :: thres_qp, tau_1_slav(3), tau_2_slav(3)
!
! ----- Project quadrature point (on master side)
!
    call projQpSl2Ma(geom, coor_qp_sl, parameters%proj_tole, &
                    coor_qp_ma, gap, tau_1_slav, tau_2_slav, norm_slav, norm_mast)
!
! ----- Evaluate shape function for displacement (slave and master)
!
    call shapeFuncDisp(geom%elem_dime, geom%nb_node_slav, geom%elem_slav_code, coor_qp_sl, &
                        shape_func_sl, dshape_func_sl)
    call shapeFuncDisp(geom%elem_dime, geom%nb_node_mast, geom%elem_mast_code, coor_qp_ma, &
                        shape_func_ma)
!
! ----- Evaluate shape function for Lagrange (slave)
!
    call shapeFuncLagr(geom%elem_dime, geom%elem_slav_code, coor_qp_sl, &
                        shape_func_lagr)
!
! ----- Evaluate Lagr_c and gamma_c at quadrature point
!
    lagr_c = evalPoly(geom%nb_lagr_c, shape_func_lagr, geom%slav_lagc_curr)
    gamma_c = evalPoly(geom%nb_lagr_c, shape_func_lagr, parameters%coef_cont) / hF
    lagrc_gap = lagr_c + gamma_c * gap
!
! ----- Contact activate at quadrature point ( H = 0 or 1 )
!
    if(parameters%type_cont == CONT_TYPE_UNIL) then
        H = Heaviside(-lagrc_gap)
        projRmVal = projRm(lagrc_gap)
    else
        H = 1.d0
        projRmVal = lagrc_gap
    end if
    l_cont_qp = (H > 0.5d0)
!
! ----- Evaluate Lagr_f and gamma_f at quadrature point
!
    vT = 0.d0
    lagr_f = 0.d0
    gamma_f = 0.d0
    projBsVal = 0.d0
    l_fric_qp = ASTER_FALSE
!
    if(parameters%l_fric) then
! ----- Define threshold for friction
        if(parameters%type_fric == FRIC_TYPE_TRES) then
            l_fric_qp = ASTER_TRUE
            thres_qp = parameters%threshold_given
        elseif(parameters%type_fric == FRIC_TYPE_NONE) then
            l_fric_qp = ASTER_FALSE
            thres_qp = 0.d0
        else
            l_fric_qp = l_cont_qp
            if(l_cont_qp) then
                if (parameters%type_fric == FRIC_TYPE_COUL) then
                    thres_qp = - parameters%threshold_given * projRmVal
                else
                    thres_qp = THRES_STICK
                end if
            else
                thres_qp = 0.d0
            end if
        end if
!
        lagr_f(1) = evalPoly(geom%nb_lagr_c, shape_func_lagr, geom%slav_lagf_curr(1,:))
        lagr_f(2) = evalPoly(geom%nb_lagr_c, shape_func_lagr, geom%slav_lagf_curr(2,:))
        gamma_f = evalPoly(geom%nb_lagr_c, shape_func_lagr, parameters%coef_fric) / hF
        lagrf_vT = lagr_f - gamma_f * vT
        projBsVal = projBs(parameters, lagrf_vT, thres_qp)
    end if
!
    if(present(dGap)) then
!
! ----- Compute d (gap(u))[v] / du
!
        call dGap_du(geom, shape_func_sl, dshape_func_sl, &
                    shape_func_ma, norm_slav, norm_mast, gap, dGap)

    end if
!
    if(present(d2Gap)) then
!
! ----- Compute d^2 (gap(u))[v, w] / du^2
!
        call d2Gap_du2(geom, shape_func_sl, dshape_func_sl, &
                        shape_func_ma, norm_slav, norm_mast, gap, d2Gap)

    end if
!
    if(present(mu_c)) then
!
! ----- Compute mu_c
!
        call testLagrC(geom, shape_func_lagr, mu_c)
    end if
!
    if(present(mu_f)) then
!
! ----- Compute mu_f
!
        call testLagrF(geom, parameters%l_fric, shape_func_lagr, mu_f)
    end if
!
    if(present(jump_t)) then
!
! ----- Compute mu_f
!
         call jump_tang(geom, shape_func_sl, shape_func_ma, tau_1_slav, tau_2_slav, jump_t)
    end if
!
!    print*, "VAL: ", lagr_c, gamma_c, gap, H, projRmVal, hF
!
end subroutine
