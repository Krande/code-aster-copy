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
subroutine niMatr(parameters, geom, matr_cont, matr_fric)
!
    use contact_module
    use contact_nitsche_module
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getInterCont.h"
#include "asterfort/getQuadCont.h"
#include "asterfort/niElemCont.h"
#include "blas/dger.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(out) :: matr_cont(MAX_NITS_DOFS, MAX_NITS_DOFS)
    real(kind=8), intent(out) :: matr_fric(MAX_NITS_DOFS, MAX_NITS_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Nitsche method) - Elementary computations
!
! Compute matrix
!
! --------------------------------------------------------------------------------------------------
!
! Out  matr_cont : matrix (only upper part)
!
! --------------------------------------------------------------------------------------------------
!
    type(ContactNitsche) :: nits
    aster_logical :: l_cont_qp, l_fric_qp
    integer(kind=8) ::  i_qp, nb_qp, total_dofs, face_dofs, slav_dofs
    real(kind=8) :: weight_sl_qp, coeff, hF
    real(kind=8) :: coor_qp_sl(2)
    real(kind=8) :: coor_qp(2, MAX_NB_QUAD), weight_qp(MAX_NB_QUAD)
    real(kind=8) :: gap, stress_nn, gamma_c, projRmVal
    real(kind=8) :: stress_t(2), vT(2), gamma_f, projBsVal(2)
    real(kind=8) :: dGap(MAX_LAGA_DOFS), d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    real(kind=8) :: dStress_nn(MAX_NITS_DOFS), dGapRenum(MAX_NITS_DOFS)
    real(kind=8) :: matr_tmp(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    integer(kind=8) :: dofsMap(54)
    integer(kind=8) :: nbPoinInte
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
!
! --------------------------------------------------------------------------------------------------
!
    matr_cont = 0.d0
    matr_fric = 0.d0
    dGapRenum = 0.d0
!
! - Mapping of dofs
    dofsMap = dofsMapping(geom)

! - Get intersection points
    call getInterCont(nbPoinInte, poinInteSlav)

! - Get quadrature (slave side)
    call getQuadCont(geom%elem_dime, &
                     geom%elem_slav_code, geom%elem_mast_code, &
                     nbPoinInte, poinInteSlav, &
                     nb_qp, coor_qp, &
                     geom%l_axis, geom%nb_node_slav, geom%coor_slav_init, &
                     weight_qp)

! - Diameter of slave side
    hF = diameter(geom%nb_node_slav, geom%coor_slav_init)

!
    call nbDofsNitsche(geom, total_dofs, face_dofs, slav_dofs)
!
! - Read material properties
!
    call getMaterialProperties(nits)
!
! - Eval stress at face nodes
!
    call evalStressNodes(geom, nits)
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
        call niElemCont(parameters, geom, nits, coor_qp_sl, hF, &
                        stress_nn, gap, gamma_c, projRmVal, l_cont_qp, &
                        stress_t, vT, gamma_f, projBsVal, l_fric_qp, &
                        dGap=dGap, d2Gap=d2Gap, dStress_nn=dStress_nn)
!
! ------ CONTACT PART (always computed)
!
        if (l_cont_qp) then
!
! ------ Compute displacement / displacement (slave and master side)
!        term: (gamma_c*H*D(gap(u))[v], D(gap(u))[du])
!
            matr_tmp = 0.d0
            coeff = weight_sl_qp*gamma_c
            b_lda = to_blas_int(MAX_LAGA_DOFS)
            b_m = to_blas_int(face_dofs)
            b_n = to_blas_int(face_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dger(b_m, b_n, coeff, dGap, b_incx, &
                      dGap, b_incy, matr_tmp, b_lda)
!
! ------ Compute displacement / displacement (slave and master side)
!        term: (H*[stress_nn + gamma_c * gap(u)]_R-, D2(gap(u))[v, du])
!
            coeff = weight_sl_qp*projRmVal
            matr_tmp(1:face_dofs, 1:face_dofs) = matr_tmp(1:face_dofs, 1:face_dofs)+coeff*d2Gap(1&
                                                 &:face_dofs, 1:face_dofs)
!
! ------ Renumbering
!
            call remappingMatr(geom, dofsMap, matr_tmp, matr_cont, 1.d0)
!
! ------ Compute displacement / displacement
!        term: (H * D(gap(u))[v], D(stress_nn)[du])
!
            call remappingVect(geom, dofsMap, dGap, dGapRenum, 1.d0)
            coeff = weight_sl_qp
            b_lda = to_blas_int(MAX_NITS_DOFS)
            b_m = to_blas_int(total_dofs)
            b_n = to_blas_int(slav_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dger(b_m, b_n, coeff, dGapRenum, b_incx, &
                      dStress_nn, b_incy, matr_cont, b_lda)
!
            if (parameters%vari_cont .ne. CONT_VARI_RAPI) then
!
! ------ Compute displacement / displacement
!        term: \theta * (H * D(stress_nn)[v], D(gap(u))[du])
!
                coeff = weight_sl_qp*parameters%vari_cont_coef
                b_lda = to_blas_int(MAX_NITS_DOFS)
                b_m = to_blas_int(slav_dofs)
                b_n = to_blas_int(total_dofs)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dger(b_m, b_n, coeff, dStress_nn, b_incx, &
                          dGapRenum, b_incy, matr_cont, b_lda)
            end if

        else
            if (parameters%vari_cont .ne. CONT_VARI_RAPI) then
!
! ------ Compute displacement / displacement
!        term: (theta * (H-1) / gamma_c * D(stress_nn)[v], D(stress_nn)[du])
!
                coeff = weight_sl_qp*parameters%vari_cont_coef*gamma_c
                b_lda = to_blas_int(MAX_NITS_DOFS)
                b_m = to_blas_int(slav_dofs)
                b_n = to_blas_int(slav_dofs)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dger(b_m, b_n, coeff, dStress_nn, b_incx, &
                          dStress_nn, b_incy, matr_cont, b_lda)
            end if
        end if
!
! ------ FRICTION PART (computed only if friction)
!
        if (parameters%l_fric) then
            if (l_fric_qp) then
                ASSERT(ASTER_FALSE)
            end if
        end if
    end do
!
end subroutine
