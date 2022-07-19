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
subroutine laVect(parameters, geom, vect)
!
use contact_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getQuadCont.h"
#include "blas/daxpy.h"
#include "contact_module.h"
#include "asterfort/laElemCont.h"
!
type(ContactParameters), intent(in) :: parameters
type(ContactGeom), intent(in) :: geom
real(kind=8), intent(inout) :: vect(MAX_CONT_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute vector
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
! IO  vect             : vector
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_qp
    integer ::  i_qp, nb_qp
    real(kind=8) :: weight_sl_qp, coeff, hF
    real(kind=8) :: coor_qp_sl(2)
    real(kind=8) :: coor_qp(2, 48), weight_qp(48)
    real(kind=8) :: gap, lagr_c, gamma_c, projRmVal
    real(kind=8) :: dGap(MAX_CONT_DOFS), mu_c(MAX_CONT_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
    vect = 0.d0
!
! - Slave node is not paired -> Special treatment
!
    if(geom%elem_slav_code == "PO1") then
        if(geom%elem_mast_code == "LAGR") then
            vect(geom%elem_dime+1) = -geom%slav_lagc_curr(1)
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
                    dGap=dGap, mu_c=mu_c)
!
! ------ Compute displacement (slave and master side)
!        term: (H*[lagr_c + gamma_c * gap(u)]_R-, D(gap(u))[v])
!
        if(l_cont_qp) then
            coeff = weight_sl_qp * projRmVal
            call daxpy(geom%nb_dofs, coeff, dGap, 1, vect, 1)
        end if
!
! ------ Compute Lagrange (slave side)
!        term: (([lagr_c + gamma_c * gap(u)]_R- - lagr_c) / gamma_c, mu_c)
!
        coeff = weight_sl_qp * (projRmVal - lagr_c) / gamma_c
        call daxpy(geom%nb_dofs, coeff, mu_c, 1, vect, 1)
    end do
!
999 continue
!
end subroutine
