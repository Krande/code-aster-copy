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
subroutine laElemCont(elem_dime, coor_qp_sl, proj_tole, &
                    nb_node_slav, elem_slav_code, slav_coor_curr,&
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    nb_lagr_c, lagc_curr, indi_lagc, gamma_c_nodes, hF, &
                    lagr_c, gap, gamma_c, projRmVal, l_cont_qp,&
                    dGap, d2Gap, mu_c)
!
use contact_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "contact_module.h"
!
integer, intent(in) :: elem_dime
integer, intent(in) :: nb_lagr_c, indi_lagc(9)
character(len=8), intent(in) :: elem_slav_code, elem_mast_code
integer, intent(in) :: nb_node_slav, nb_node_mast
real(kind=8), intent(in) :: slav_coor_curr(3, 9), mast_coor_curr(3, 9)
real(kind=8), intent(in) :: coor_qp_sl(2), hF
real(kind=8), intent(in) :: proj_tole, gamma_c_nodes(4), lagc_curr(4)
real(kind=8), intent(out) :: lagr_c, gap, gamma_c, projRmVal
aster_logical, intent(out) :: l_cont_qp
real(kind=8), intent(out), optional :: dGap(MAX_CONT_DOFS)
real(kind=8), intent(out), optional :: d2Gap(MAX_CONT_DOFS, MAX_CONT_DOFS)
real(kind=8), intent(out), optional :: mu_c(MAX_CONT_DOFS)
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
    real(kind=8) :: norm_slav(3), norm_mast(3), H, coor_qp_ma(2)
!
! ----- Project quadrature point (on master side)
!
    call projQpSl2Ma(elem_dime, coor_qp_sl, proj_tole, &
                    nb_node_slav, elem_slav_code, slav_coor_curr,&
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    coor_qp_ma, gap, norm_slav, norm_mast)
!
! ----- Evaluate shape function for displacement (slave and master)
!
    call shapeFuncDisp(elem_dime, nb_node_slav, elem_slav_code, coor_qp_sl, &
                        shape_func_sl, dshape_func_sl)
    call shapeFuncDisp(elem_dime, nb_node_mast, elem_mast_code, coor_qp_ma, &
                        shape_func_ma)
!
! ----- Evaluate shape function for Lagrange (slave)
!
    call shapeFuncLagr(elem_dime, elem_slav_code, coor_qp_sl, &
                        shape_func_lagr)
!
! ----- Evaluate Lagr_c and gamma_c at quadrature point
!
    lagr_c = evalPoly(nb_lagr_c, shape_func_lagr, lagc_curr)
    gamma_c = evalPoly(nb_lagr_c, shape_func_lagr, gamma_c_nodes) / hF
    projRmVal = projRm(lagr_c + gamma_c * gap)
!
! ----- Contact activate at quadrature point ( H = 0 or 1 )
!
    H = Heaviside(-(lagr_c + gamma_c * gap))
    l_cont_qp = (H > 0.5d0)
!
    if(present(dGap)) then
!
! ----- Compute d (gap(u))[v] / du
!
        call dGap_du(elem_dime, nb_node_slav, shape_func_sl, dshape_func_sl, &
                    nb_node_mast, shape_func_ma, &
                    indi_lagc, norm_slav, norm_mast, gap, dGap)

    end if
!
    if(present(d2Gap)) then
!
! ----- Compute d^2 (gap(u))[v, w] / du^2
!
        call d2Gap_du2(elem_dime, nb_node_slav, shape_func_sl, nb_node_mast, shape_func_ma, &
            indi_lagc, norm_slav, norm_mast, d2Gap)

    end if
!
    if(present(mu_c)) then
!
! ----- Compute mu_c
!
        call testLagrC(elem_dime, nb_lagr_c, shape_func_lagr, indi_lagc, mu_c)
    end if
!
!    print*, "VAL: ", lagr_c, gamma_c, gap, H, projRmVal, hF
!
end subroutine
