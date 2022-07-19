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
subroutine laVect(elem_dime   , l_axis        , nb_dofs, &
                  nb_lagr_c   , indi_lagc     , lagc_curr, &
                  gamma_c_nodes, &
                  nb_node_slav, elem_slav_code, slav_coor_init, slav_coor_curr,&
                  nb_node_mast, elem_mast_code, mast_coor_curr,&
                  proj_tole, vect)
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
integer, intent(in) :: elem_dime
aster_logical, intent(in) :: l_axis
integer, intent(in) :: nb_lagr_c, indi_lagc(9), nb_dofs
character(len=8), intent(in) :: elem_slav_code, elem_mast_code
integer, intent(in) :: nb_node_slav, nb_node_mast
real(kind=8), intent(in) :: slav_coor_curr(3, 9), slav_coor_init(3,9)
real(kind=8), intent(in) :: mast_coor_curr(3, 9)
real(kind=8), intent(in) :: proj_tole, gamma_c_nodes(4), lagc_curr(4)
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
    if(elem_slav_code == "PO1") then
        if(elem_mast_code == "LAGR") then
            vect(elem_dime+1) = -lagc_curr(1)
        else
            if(elem_mast_code .ne. "NOLAGR") then
                ASSERT(ASTER_FALSE)
            end if
        end if
!
        go to 999
    end if
!
! - Get quadrature (slave side)
!
    call getQuadCont(elem_dime, l_axis, nb_node_slav, elem_slav_code, slav_coor_init,&
                     elem_mast_code, nb_qp, coor_qp, weight_qp )
!
! - Diameter of slave side
!
    hF = diameter(nb_node_slav, slav_coor_init)
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
        call laElemCont(elem_dime, coor_qp_sl, proj_tole, &
                    nb_node_slav, elem_slav_code, slav_coor_curr,&
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    nb_lagr_c, lagc_curr, indi_lagc, gamma_c_nodes, hF, &
                    lagr_c, gap, gamma_c, projRmVal, l_cont_qp,&
                    dGap=dGap, mu_c=mu_c)
!
! ------ Compute displacement (slave and master side)
!        term: (H*[lagr_c + gamma_c * gap(u)]_R-, D(gap(u))[v])
!
        if(l_cont_qp) then
            coeff = weight_sl_qp * projRmVal
            call daxpy(nb_dofs, coeff, dGap, 1, vect, 1)
        end if
!
! ------ Compute Lagrange (slave side)
!        term: (([lagr_c + gamma_c * gap(u)]_R- - lagr_c) / gamma_c, mu_c)
!
        coeff = weight_sl_qp * (projRmVal - lagr_c) / gamma_c
        call daxpy(nb_dofs, coeff, mu_c, 1, vect, 1)
    end do
!
999 continue
!
end subroutine
