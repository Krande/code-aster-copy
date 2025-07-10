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
subroutine peMatr_ct_pr(parameters, geom, matr_cont, matr_fric)
!
    use contact_type
    use contact_module
    use contact_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/getInterCont.h"
#include "asterfort/getQuadCont.h"
#include "asterfort/laElemCont.h"
#include "blas/dgemm.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: matr_cont(MAX_PENA_DOFS, MAX_PENA_DOFS)
    real(kind=8), intent(inout) :: matr_fric(MAX_PENA_DOFS, MAX_PENA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Penalization method) - Elementary computations
!
! Compute matrices
!
! WARNING
! Note that we build the transpose of the Jacobian - it allows the computer expressions to be
! closer to the mathematical ones. The matrix is then transposed at the end of the subroutine.
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : numerical parameters
! In  geom             : geometrical information
! IO  matr_cont        : matrix (only upper part)
! IO  matr_fric        : matrix
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_qp, l_fric_qp, l_print
    integer(kind=8) :: i_qp, nb_qp, i_dof, i_zeta, i_cmp, j_dof, nbPoinInte
    real(kind=8) :: weight_sl_qp, coeff, hF
    real(kind=8) :: coor_qp_sl(2), tau_slav(3, 2)
    real(kind=8) :: coor_qp(2, MAX_NB_QUAD), weight_qp(MAX_NB_QUAD)
    real(kind=8) :: dts_ns(MAX_PENA_DOFS, 2)
    real(kind=8) :: gap, gamma_c, projRmVal, dNs(MAX_PENA_DOFS, 3)
    real(kind=8) :: gamma_f, dGap(MAX_PENA_DOFS)
    real(kind=8) :: norm_slav(3), jump_v(MAX_PENA_DOFS, 3)
    real(kind=8) :: mu_c_ns(MAX_PENA_DOFS, 3)
    real(kind=8) :: speed(3), H, dfunc_ma(3, MAX_PENA_DOFS, 2), dZetaM(MAX_PENA_DOFS, 2)
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    real(kind=8), pointer :: dfunc_dzeta(:, :, :) => null()
    blas_int :: b_dime, b_nb_dofs
    blas_int :: b_1, b_3, b_MAX_PENA_DOFS
!
! --------------------------------------------------------------------------------------------------
!
    b_MAX_PENA_DOFS = to_blas_int(MAX_PENA_DOFS)
    b_1 = 1
    b_3 = 3
    l_print = ASTER_FALSE
!
! - Large arrays allocated on the heap rather than on the stack
    allocate (dfunc_dzeta(3, MAX_PENA_DOFS, MAX_PENA_DOFS))
!
    matr_cont = 0.d0
    matr_fric = 0.d0
!
! - Slave node is not paired -> Special treatment
!
    if (geom%elem_slav_code == "PO1") then
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
                        gap, gamma_c, projRmVal, l_cont_qp, &
                        gamma_f, l_fric_qp, norm_slav=norm_slav, dGap=dGap, &
                        jump_v=jump_v, dNs=dNs, &
                        tau_slav=tau_slav, dts_ns=dts_ns, &
                        speed=speed, dfunc_ma=dfunc_ma, &
                        dZetaM=dZetaM)
!
! - Contact indicator
!
        H = 0.d0
        if (l_cont_qp) H = 1.d0
!
! ------ CONTACT PART (always computed)
!
! ------ Compute displacement / Lagrange and Lagrange / displacement
!        term: -(dlagr_c n^s, -v^m + v^s)
!
        ! write (6, *) 'l_cont_qp', l_cont_qp
        if (l_cont_qp) then
            coeff = -weight_sl_qp*gamma_c
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            b_dime = to_blas_int(geom%elem_dime)
            mu_c_ns = 0.d0
            do i_dof = 1, geom%nb_dofs
                mu_c_ns(i_dof, :) = dGap(i_dof)*norm_slav
            end do
            call dgemm('N', 'T', b_nb_dofs, b_nb_dofs, b_dime, coeff, &
                       mu_c_ns, b_MAX_PENA_DOFS, jump_v, b_MAX_PENA_DOFS, 1.d0, &
                       matr_cont, b_MAX_PENA_DOFS)
            if (l_print) then
                write (6, *) '------------------------------------------------------------'
                do i_dof = 1, b_nb_dofs
                    write (6, 10) 'matr_cont0=', matr_cont(1:b_nb_dofs, i_dof)
                end do
            end if
!
!        term: -([gamma_c*gap]_ Dn^s[du], -v^m + v^s)
!
            coeff = -weight_sl_qp*projRmVal
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            b_dime = to_blas_int(geom%elem_dime)
            call dgemm('N', 'T', b_nb_dofs, b_nb_dofs, b_dime, coeff, &
                       dNs, b_MAX_PENA_DOFS, jump_v, b_MAX_PENA_DOFS, 1.d0, &
                       matr_cont, b_MAX_PENA_DOFS)
            if (l_print) then
                write (6, *) '------------------------------------------------------------'
                do i_dof = 1, b_nb_dofs
                    write (6, 10) 'matr_cont1=', matr_cont(1:b_nb_dofs, i_dof)
                end do
            end if
!
!        term: ([gamma_c*gap]_ n^s, Dv^m[du])  (master side)
!
!           Dv^m[du] = (d v^m/d zeta) * D zeta^m[du]
            dfunc_dzeta = 0.d0
            do i_cmp = 1, geom%elem_dime
                do i_dof = 1, geom%nb_dofs
                    do j_dof = 1, geom%nb_dofs
                        do i_zeta = 1, geom%elem_dime-1
                            dfunc_dzeta(i_cmp, i_dof, j_dof) = dfunc_dzeta(i_cmp, i_dof, j_dof)+ &
                                                               dfunc_ma(i_cmp, j_dof, i_zeta)* &
                                                               dZetaM(i_dof, i_zeta)
                        end do
                    end do
                end do
            end do
!           ([gamma_c*gap]_*norm_slav, Dv^m[du])
            coeff = weight_sl_qp
            do i_cmp = 1, geom%elem_dime
                matr_cont = matr_cont+coeff*projRmVal*norm_slav(i_cmp)*dfunc_dzeta(i_cmp, :, :)
            end do
            if (l_print) then
                write (6, *) '------------------------------------------------------------'
                do i_dof = 1, b_nb_dofs
                    write (6, 10) 'matr_cont2=', matr_cont(1:b_nb_dofs, i_dof)
                end do
            end if
!
!
        end if
!
!
    end do
!
999 continue
10  FORMAT(a, *(F6.3))
!
! - deallocate large array
    deallocate (dfunc_dzeta)
!
! ------ Transpose (see warning at head of the routine)
!
    matr_cont = transpose(matr_cont)
!
end subroutine
