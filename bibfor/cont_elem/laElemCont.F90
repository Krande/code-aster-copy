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
type(ContactParameters), intent(in) :: parameters
type(ContactGeom), intent(in) :: geom
real(kind=8), intent(in) :: coor_qp_sl(2), hF
real(kind=8), intent(out) :: lagr_c, gap, gamma_c, projRmVal
aster_logical, intent(out) :: l_cont_qp
real(kind=8), intent(out), optional :: dGap(MAX_LAGA_DOFS)
real(kind=8), intent(out), optional :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
real(kind=8), intent(out), optional :: mu_c(MAX_LAGA_DOFS)
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
    real(kind=8) :: norm_slav(3), norm_mast(3), H, coor_qp_ma(2), lagr_gap
!
! ----- Project quadrature point (on master side)
!
    call projQpSl2Ma(geom, coor_qp_sl, parameters%proj_tole, &
                    coor_qp_ma, gap, norm_slav, norm_mast)
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
    lagr_gap = lagr_c + gamma_c * gap
!
! ----- Contact activate at quadrature point ( H = 0 or 1 )
!
    if(parameters%type_cont == CONT_TYPE_UNIL) then
        H = Heaviside(-lagr_gap)
        projRmVal = projRm(lagr_gap)
    else
        H = 1.d0
        projRmVal = lagr_gap
    end if
    l_cont_qp = (H > 0.5d0)
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
        call d2Gap_du2(geom, shape_func_sl, shape_func_ma, &
                        norm_slav, norm_mast, d2Gap)

    end if
!
    if(present(mu_c)) then
!
! ----- Compute mu_c
!
        call testLagrC(geom, shape_func_lagr, mu_c)
    end if
!
!    print*, "VAL: ", lagr_c, gamma_c, gap, H, projRmVal, hF
!
end subroutine
