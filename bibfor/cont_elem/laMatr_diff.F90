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
subroutine laMatr_diff(parameters, geom, matr_cont, matr_fric)
!
    use contact_type
    use contact_module
    use contact_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/getInterCont.h"
#include "asterfort/getQuadCont.h"
#include "asterfort/laVect.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(inout) :: geom
    real(kind=8), intent(inout) :: matr_cont(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    real(kind=8), intent(inout) :: matr_fric(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute matrix by finite differences
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : parameters of the method
! In  geom             : input data
! IO  matr_cont        : matrix
! IO  matr_fric        : matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_qp, nb_qp, i_lagr, i_node, i_dim, index, order, nbPoinInte
    real(kind=8) :: hF, coor_qp_sl(2), eps, norm_slav(3)
    real(kind=8) :: coor_qp(2, MAX_NB_QUAD), weight_qp(MAX_NB_QUAD)
    real(kind=8) :: vect_cont(MAX_LAGA_DOFS), vect_fric(MAX_LAGA_DOFS)
    real(kind=8) :: vect_cont_p(MAX_LAGA_DOFS), vect_fric_p(MAX_LAGA_DOFS)
    real(kind=8) :: vect_cont_m(MAX_LAGA_DOFS), vect_fric_m(MAX_LAGA_DOFS)
    real(kind=8) :: coor_save, depl_save, pair_save, lagc_save, lagf_save
    real(kind=8) :: tau_slav(3, 2), delta
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    character(len=8), dimension(3) :: l_dof
    character(len=8) :: ref
!
! --------------------------------------------------------------------------------------------------
!
    matr_cont = 0.d0
    matr_fric = 0.d0
    order = 1
    l_dof = [character(len=8) :: "DX", "DY", "DZ"]
    ref = 'ref'
!
! - Slave node is not paired -> Special treatment
!
    if (geom%elem_slav_code == "PO1") then
        if (geom%elem_mast_code == "LAGR") then
            matr_cont(geom%elem_dime+1, geom%elem_dime+1) = 1.d0
            if (parameters%l_fric) then
                matr_fric(geom%elem_dime+2, geom%elem_dime+2) = 1.d0
                if (geom%elem_dime == 3) then
                    matr_fric(geom%elem_dime+3, geom%elem_dime+3) = 1.d0
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
! - Compute outward slave normal (to have some idea...)
!
! - Get first quadrature point (slave side)
    i_qp = 1
    coor_qp_sl(1:2) = coor_qp(1:2, i_qp)
    call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_pair, &
                coor_qp_sl(1), coor_qp_sl(2), norm_slav, tau_slav(1:3, 1), tau_slav(1:3, 2))

!
! --- Compute reference contact residual
!
    call laVect(parameters, geom, vect_cont, vect_fric)
!
! - Loop over dofs
!
    eps = 1.d-6*hF
!
! --- Slave side
!
    i_lagr = 1
    index = 1
    do i_node = 1, geom%nb_node_slav
        ! Displacement
        do i_dim = 1, geom%elem_dime
            ! Save vec
            coor_save = geom%coor_slav_curr(i_dim, i_node)
            depl_save = geom%depl_slav_curr(i_dim, i_node)
            pair_save = geom%coor_slav_pair(i_dim, i_node)
! ------------------------Order 1--------------------------------------
            if (order .eq. 1) then
                ! Some tricks to perturb displ to activate contact
                if (abs(norm_slav(i_dim)) .gt. eps) then
                    geom%coor_slav_curr(i_dim, i_node) = coor_save+eps*norm_slav(i_dim)
                    geom%depl_slav_curr(i_dim, i_node) = depl_save+eps*norm_slav(i_dim)
                    geom%coor_slav_pair(i_dim, i_node) = pair_save+eps*norm_slav(i_dim)
                else
                    geom%coor_slav_curr(i_dim, i_node) = coor_save+eps
                    geom%depl_slav_curr(i_dim, i_node) = depl_save+eps
                    geom%coor_slav_pair(i_dim, i_node) = pair_save+eps
                end if
                ! Compute perturbed residuals
                call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                ! Compute matrices
                delta = geom%depl_slav_curr(i_dim, i_node)-depl_save
                ! write (6, *) '*delta*', delta
                matr_cont(:, index) = (vect_cont_p-vect_cont)/delta
                matr_fric(:, index) = (vect_fric_p-vect_fric)/delta
            else
! ------------------------Order 2--------------------------------------
                ! Perturb displ +
                if (abs(norm_slav(i_dim)) .gt. eps) then
                    geom%coor_slav_curr(i_dim, i_node) = coor_save-eps*norm_slav(i_dim)
                    geom%depl_slav_curr(i_dim, i_node) = depl_save-eps*norm_slav(i_dim)
                    geom%coor_slav_pair(i_dim, i_node) = pair_save-eps*norm_slav(i_dim)
                else
                    geom%coor_slav_curr(i_dim, i_node) = coor_save+eps
                    geom%depl_slav_curr(i_dim, i_node) = depl_save+eps
                    geom%coor_slav_pair(i_dim, i_node) = pair_save+eps
                end if
                ! Compute perturbed residuals
                call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                ! Perturb displ
                if (abs(norm_slav(i_dim)) .gt. eps) then
                    geom%coor_slav_curr(i_dim, i_node) = coor_save+eps*norm_slav(i_dim)
                    geom%depl_slav_curr(i_dim, i_node) = depl_save+eps*norm_slav(i_dim)
                    geom%coor_slav_pair(i_dim, i_node) = pair_save+eps*norm_slav(i_dim)
                else
                    geom%coor_slav_curr(i_dim, i_node) = coor_save-eps
                    geom%depl_slav_curr(i_dim, i_node) = depl_save-eps
                    geom%coor_slav_pair(i_dim, i_node) = pair_save-eps
                end if
                ! Compute perturbed residuals
                call laVect(parameters, geom, vect_cont_m, vect_fric_m)
                ! Compute matrices
                delta = coor_save-geom%coor_slav_curr(i_dim, i_node)
                matr_cont(:, index) = (vect_cont_p-vect_cont_m)/delta/2
                matr_fric(:, index) = (vect_fric_p-vect_fric_m)/delta/2
            end if
! ----------------------------------------------------------------------
            ! Restore vec
            geom%coor_slav_curr(i_dim, i_node) = coor_save
            geom%depl_slav_curr(i_dim, i_node) = depl_save
            geom%coor_slav_pair(i_dim, i_node) = pair_save
            index = index+1
        end do
! ==================================================================================================
        ! Does this node handle Lagrange multipliers ?
        if (geom%indi_lagc(i_node) > 0) then
            ! Contact Lagrange multiplier
            ! Save lagr_c
            lagc_save = geom%lagc_slav_curr(i_lagr)
! ------------------------Order 1--------------------------------------
            if (order .eq. 1) then
                ! Perturb lagr_c
                if (abs(lagc_save) .gt. eps) then
                    geom%lagc_slav_curr(i_lagr) = lagc_save*(1.d0-eps)
                else
                    geom%lagc_slav_curr(i_lagr) = -eps
                end if
                ! Compute perturbed residuals
                call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                ! Compute matrices
                delta = geom%lagc_slav_curr(i_lagr)-lagc_save
                matr_cont(:, index) = (vect_cont_p-vect_cont)/delta
                matr_fric(:, index) = (vect_fric_p-vect_fric)/delta
            else
! ------------------------Order 2--------------------------------------
                ! Perturb lagr_c+
                if (abs(lagc_save) .gt. eps) then
                    geom%lagc_slav_curr(i_lagr) = lagc_save*(1.d0-eps)
                else
                    geom%lagc_slav_curr(i_lagr) = -eps
                end if
                ! Compute perturbed residuals
                call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                ! Perturb lagr_c-
                if (abs(lagc_save) .gt. eps) then
                    geom%lagc_slav_curr(i_lagr) = lagc_save*(1.d0+eps)
                else
                    geom%lagc_slav_curr(i_lagr) = +eps
                end if
                ! Compute perturbed residuals
                call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                ! Compute matrices
                delta = lagc_save-geom%lagc_slav_curr(i_lagr)
                matr_cont(:, index) = (vect_cont_p-vect_cont_m)/delta/2
                matr_fric(:, index) = (vect_fric_p-vect_fric_m)/delta/2
            end if
! ---------------------------------------------------------------------
            ! Restore lagr_c
            geom%lagc_slav_curr(i_lagr) = lagc_save
            index = index+1
            if (parameters%l_fric) then
                ! Friction Lagrange multipliers
                do i_dim = 1, geom%elem_dime-1
                    ! Save lagr_f
                    lagf_save = geom%lagf_slav_curr(i_dim, i_lagr)
! ------------------------Order 1--------------------------------------
                    if (order .eq. 1) then
                        ! Perturb lagr_f
                        if (abs(lagf_save) .gt. eps) then
                            geom%lagf_slav_curr(i_dim, i_lagr) = lagf_save*(1.d0-eps)
                            ! write (6, *) 'lagf != 0', lagf_save
                        else
                            geom%lagf_slav_curr(i_dim, i_lagr) = -eps
                            ! write (6, *) 'lagf == 0', lagf_save
                        end if
                        ! Compute perturbed residuals
                        call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                        ! Compute matrices
                        delta = geom%lagf_slav_curr(i_dim, i_lagr)-lagf_save
                        matr_cont(:, index) = (vect_cont_p-vect_cont)/delta
                        matr_fric(:, index) = (vect_fric_p-vect_fric)/delta
                    else
! ------------------------Order 2--------------------------------------
                        ! Perturb lagr_f+
                        if (abs(lagf_save) .gt. eps) then
                            geom%lagf_slav_curr(i_dim, i_lagr) = lagf_save*(1.d0-eps)
                        else
                            geom%lagf_slav_curr(i_dim, i_lagr) = -eps
                        end if
                        ! Compute perturbed residuals
                        call laVect(parameters, geom, vect_cont_p, vect_fric_p)
                        ! Perturb lagr_f-
                        if (abs(lagf_save) .gt. eps) then
                            geom%lagf_slav_curr(i_dim, i_lagr) = lagf_save*(1.d0+eps)
                        else
                            geom%lagf_slav_curr(i_dim, i_lagr) = +eps
                        end if
                        ! Compute perturbed residuals
                        call laVect(parameters, geom, vect_cont_m, vect_fric_m)
                        ! Compute matrices
                        delta = lagf_save-geom%lagf_slav_curr(i_dim, i_lagr)
                        matr_cont(:, index) = (vect_cont_p-vect_cont_m)/delta/2
                        matr_fric(:, index) = (vect_fric_p-vect_fric_m)/delta/2
                    end if
! --------------------------------------------------------------------
                    ! Restore lagr_f
                    geom%lagf_slav_curr(i_dim, i_lagr) = lagf_save
                    index = index+1
                end do
            end if
            i_lagr = i_lagr+1
        end if
    end do
!
! --- Master side
!
    do i_node = 1, geom%nb_node_mast
        ! Displacement
        do i_dim = 1, geom%elem_dime
            ! Save vec
            coor_save = geom%coor_mast_curr(i_dim, i_node)
            depl_save = geom%depl_mast_curr(i_dim, i_node)
            pair_save = geom%coor_mast_pair(i_dim, i_node)
            ! Some tricks to perturb displ to activate contact
            if (abs(norm_slav(i_dim)) .gt. eps) then
                geom%coor_mast_curr(i_dim, i_node) = coor_save-eps*norm_slav(i_dim)
                geom%depl_mast_curr(i_dim, i_node) = depl_save-eps*norm_slav(i_dim)
                geom%coor_mast_pair(i_dim, i_node) = pair_save-eps*norm_slav(i_dim)
            else
                geom%coor_mast_curr(i_dim, i_node) = coor_save-eps
                geom%depl_mast_curr(i_dim, i_node) = depl_save-eps
                geom%coor_mast_pair(i_dim, i_node) = pair_save-eps
            end if
!
            ! Compute perturbed residuals
            call laVect(parameters, geom, vect_cont_p, vect_fric_p)
            ! Compute matrices
            delta = geom%coor_mast_curr(i_dim, i_node)-coor_save
            matr_cont(:, index) = (vect_cont_p-vect_cont)/delta
            matr_fric(:, index) = (vect_fric_p-vect_fric)/delta
            ! Restore vec
            geom%coor_mast_curr(i_dim, i_node) = coor_save
            geom%depl_mast_curr(i_dim, i_node) = depl_save
            geom%coor_mast_pair(i_dim, i_node) = pair_save
            index = index+1
        end do
    end do
!
!
999 continue
!
end subroutine
