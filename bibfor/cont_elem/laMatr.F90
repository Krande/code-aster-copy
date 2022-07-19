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
subroutine laMatr(parameters, geom, matr)
!
use contact_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getQuadCont.h"
#include "blas/dger.h"
#include "blas/dsyr.h"
#include "contact_module.h"
#include "asterfort/laElemCont.h"
!
type(ContactParameters), intent(in) :: parameters
type(ContactGeom), intent(in) :: geom
real(kind=8), intent(inout) :: matr(MAX_CONT_DOFS, MAX_CONT_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute matrix
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
! IO  matr             : matrix (only upper part)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_qp
    integer ::  i_qp, nb_qp
    real(kind=8) :: weight_sl_qp, coeff, hF
    real(kind=8) :: coor_qp_sl(2)
    real(kind=8) :: coor_qp(2, 48), weight_qp(48)
    real(kind=8) :: gap, lagr_c, gamma_c, projRmVal
    real(kind=8) :: dGap(MAX_CONT_DOFS), mu_c(MAX_CONT_DOFS), d2Gap(MAX_CONT_DOFS, MAX_CONT_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
    matr = 0.d0
!
! - Slave node is not paired -> Special treatment
!
    if(geom%elem_slav_code == "PO1") then
        if(geom%elem_mast_code == "LAGR") then
            matr(geom%elem_dime+1, geom%elem_dime+1) = 1.d0
        else
            if(geom%elem_mast_code .ne. "NOLAGR") then
                ASSERT(ASTER_FALSE)
            end if
        end if
!
        go to 999
    end if
!
! - Get quadrature (slave side)
!
    call getQuadCont(geom%elem_dime, geom%l_axis, geom%nb_node_slav, geom%elem_slav_code, &
                     geom%slav_coor_init, geom%elem_mast_code, nb_qp, coor_qp, weight_qp )
!
! - Diameter of slave side
!
    hF = diameter(geom%nb_node_slav, geom%slav_coor_init)
!
! - Loop on quadrature points
!
    do i_qp = 1, nb_qp
!
! ----- Get current quadrature point (slave side)
!
        coor_qp_sl(1:2) = coor_qp(1:2, i_qp)
        weight_sl_qp = weight_qp(i_qp)
!
! ----- Compute contact quantities
!
        call laElemCont(parameters, geom, coor_qp_sl, hF, &
                    lagr_c, gap, gamma_c, projRmVal, l_cont_qp,&
                    dGap=dGap, d2Gap=d2Gap, mu_c=mu_c)
!
        if(l_cont_qp) then
!
! ------ Compute displacement / displacement (slave and master side)
!        term: (gamma_c*H*D(gap(u))[v], D(gap(u))[du])
!
            coeff = weight_sl_qp * gamma_c
            call dger(geom%nb_dofs, geom%nb_dofs, coeff, dGap, 1, dGap, 1, matr, MAX_CONT_DOFS)
!
! ------ Compute displacement / displacement (slave and master side)
!        term: (H*[lagr_c + gamma_c * gap(u)]_R-, D2(gap(u))[v, du]) -> not implemented
!
            coeff = weight_sl_qp * projRmVal
            matr(1:geom%nb_dofs, 1:geom%nb_dofs) = matr(1:geom%nb_dofs, 1:geom%nb_dofs) + &
                coeff * d2Gap(1:geom%nb_dofs, 1:geom%nb_dofs)
!
! ------ Compute displacement / Lagrange and Lagrange / displacement
!        term: (H * D(gap(u))[v], dlagr_c) -> Upper part
!        term: (H * mu_c,  D(gap(u))[du]) -> Lower part
!
            coeff = weight_sl_qp
            call dger(geom%nb_dofs, geom%nb_dofs, coeff, dGap, 1, mu_c, 1, matr, MAX_CONT_DOFS)
            call dger(geom%nb_dofs, geom%nb_dofs, coeff, mu_c, 1, dGap, 1, matr, MAX_CONT_DOFS)
        else
!
! ------ Compute Lagrange / Lagrange (slave side)
!        term: ((H-1) / gamma_c * mu_c, dlagr_c) = (-1/ gamma_c * mu_c, dlagr_c) since H = 0
!
            coeff = -weight_sl_qp / gamma_c
            call dger(geom%nb_dofs, geom%nb_dofs, coeff, mu_c, 1, mu_c, 1, matr, MAX_CONT_DOFS)
!
        end if
    end do
!
999 continue
!
end subroutine
