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
!
subroutine laVect_cf_std(parameters, geom, vect_cont, vect_fric, k_diff)
!
    use contact_module
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getInterCont.h"
#include "asterfort/getQuadCont.h"
#include "asterfort/laElemCont.h"
#include "blas/daxpy.h"
#include "blas/dgemv.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: vect_cont(MAX_LAGA_DOFS), vect_fric(MAX_LAGA_DOFS)
    character(len=8), intent(in), optional :: k_diff
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute vector
!
! -----
!
! The implementation follows Poulios' & Renard's paper : An unconstrained integral approximation
! of large sliding frictional contact between deformable solids. Computers & Structures, 2015
! The Lagrange multipliers are expressed in global coordinates
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : numerical parameters
! In  geom             : geometrical information
! IO  vect_cont        : vector for contact
! IO  vect_fric        : vector for friction
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_qp, l_fric_qp
    integer(kind=8) :: i_qp, nb_qp, nbPoinInte
    real(kind=8) :: weight_sl_qp, coeff, hF, lagr_f3(3), mu_f3(MAX_LAGA_DOFS, 3)
    real(kind=8) :: coor_qp_sl(2), metricTens(2, 2)
    real(kind=8) :: coor_qp(2, MAX_NB_QUAD), weight_qp(MAX_NB_QUAD)
    real(kind=8) :: gap, lagr_c, gamma_c, projRmVal, term_f3(3)
    real(kind=8) :: lagr_f(2), gamma_f, projBsVal(3), projBsVal2(2), term_f(2)
    real(kind=8) :: dGap(MAX_LAGA_DOFS), mu_c(MAX_LAGA_DOFS)
    real(kind=8) :: mu_f(MAX_LAGA_DOFS, 2), dvT_du(MAX_LAGA_DOFS, 3)
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    blas_int :: b_1, b_nb_dofs, b_dime_m1, b_dime
    blas_int :: b_MAX_LAGA_DOFS
!
! --------------------------------------------------------------------------------------------------
!
    b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
    b_1 = 1
!
    vect_cont = 0.d0
    vect_fric = 0.d0
!
    if (present(k_diff)) then
        ! to avoid compilation warning of dummy argument
    end if
!
! - Slave node is not paired -> Special treatment
!
    if (geom%elem_slav_code == "PO1") then
        if (geom%elem_mast_code == "LAGR") then
            vect_cont(geom%elem_dime+1) = -geom%lagc_slav_curr(1)
            if (parameters%l_fric) then
                vect_fric(geom%elem_dime+2) = -geom%lagf_slav_curr(1, 1)
                if (geom%elem_dime == 3) then
                    vect_fric(geom%elem_dime+3) = -geom%lagf_slav_curr(2, 1)
                end if
            end if
        else
            if (geom%elem_mast_code .ne. "NOLAGR") then
                ASSERT(ASTER_FALSE)
            end if
        end if
!
        go to 999
    end if
!
! - Get intersection points
!
    call getInterCont(nbPoinInte, poinInteSlav)
!
! - Get quadrature (slave side)
!
    call getQuadCont(geom%elem_dime, &
                     geom%elem_slav_code, geom%elem_mast_code, &
                     nbPoinInte, poinInteSlav, &
                     nb_qp, coor_qp, &
                     geom%l_axis, geom%nb_node_slav, geom%coor_slav_init, &
                     weight_qp)

!
! - Diameter of slave side
!
    hF = diameter(geom%nb_node_slav, geom%coor_slav_init)
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
                        gap, gamma_c, projRmVal, l_cont_qp, gamma_f, l_fric_qp, &
                        lagr_c=lagr_c, lagr_f=lagr_f, dGap=dGap, &
                        mu_c=mu_c, mu_f=mu_f, dvT=dvT_du, projBsVal=projBsVal, &
                        projBsVal2=projBsVal2, metricTens=metricTens, mu_f3=mu_f3, &
                        lagr_f3=lagr_f3)
!
! ------ CONTACT PART (always computed)
!
        if (l_cont_qp) then
!
! ------ Compute displacement (slave and master side)
!        term: (H*[lagr_c + gamma_c * gap(u)]_R-, D(gap(u))[v])
!
            coeff = weight_sl_qp*projRmVal
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            call daxpy(b_nb_dofs, coeff, dGap, b_1, vect_cont, &
                       b_1)
        end if
!
! ------ Compute Lagrange (slave side)
!        term: (([lagr_c + gamma_c * gap(u)]_R- - lagr_c) / gamma_c, mu_c)
!
        coeff = weight_sl_qp*(projRmVal-lagr_c)/gamma_c
        b_nb_dofs = to_blas_int(geom%nb_dofs)
        call daxpy(b_nb_dofs, coeff, mu_c, b_1, vect_cont, &
                   b_1)
!
!
! ------ FRICTION PART (computed only if friction)
!
        if (parameters%l_fric) then
!
! ------ Compute displacement (slave and master side)
!        term: -([lagr_f - gamma_f * vT(u)]_Bs, D(vT(u))[v])
!
            if (l_fric_qp) then
                coeff = -weight_sl_qp
                b_nb_dofs = to_blas_int(geom%nb_dofs)
                b_dime = to_blas_int(geom%elem_dime)
                call dgemv('N', b_nb_dofs, b_dime, coeff, dvT_du, &
                           b_MAX_LAGA_DOFS, projBsVal, b_1, 1.d0, vect_fric, &
                           b_1)
            end if
!
! ------ Compute Lagrange (slave side)
!        term: (([lagr_f - gamma_f * vT(u)]_Bs - lagr_f) / gamma_f, mu_f)
!
            coeff = weight_sl_qp
            term_f = matmul(metricTens, (projBsVal2-lagr_f)/gamma_f)
            term_f3 = (projBsVal-lagr_f3)/gamma_f
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            b_dime_m1 = to_blas_int(geom%elem_dime-1)
            call dgemv('N', b_nb_dofs, b_dime_m1, coeff, mu_f, &
                       b_MAX_LAGA_DOFS, term_f, b_1, 1.d0, vect_fric, &
                       b_1)
        end if
    end do
!
999 continue
!
!
end subroutine
