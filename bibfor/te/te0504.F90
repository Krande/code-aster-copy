! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine te0504(option, nomte)
!
    use HHO_basis_module
    use HHO_eval_module
    use HHO_init_module
    use HHO_matrix_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_SmallStrainMeca_module
    use HHO_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/rcangm.h"
#include "asterfort/readVector.h"
#include "asterfort/sigtmc.h"
#include "asterfort/writeVector.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: HHO
!
! Option: SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_Cell) :: hhoCell
    type(HHO_Data) :: hhoData
    type(HHO_matrix) :: gradsym
    type(HHO_Quadrature) :: hhoQuadCellRigi
!
    integer(kind=8) :: cbs, fbs, total_dofs, npg, ipg, gbs, gbs_sym, cbs_cmp
    integer(kind=8) :: nbsig, jvMater, faces_dofs, gbs_cmp, i, j
    character(len=4), parameter :: fami = 'RIGI'
    real(kind=8) :: time, sigma(6), weight, coorpg(3), anglNaut(3)
    real(kind=8) :: E_coeff(MSIZE_CELL_MAT), Eps(6)
    real(kind=8) :: dmat(6, 6), BSCEval(MSIZE_CELL_SCAL)
    real(kind=8) :: sigmVarc(6*MAX_QP_CELL), sief(6*MAX_QP_CELL)
    real(kind=8) :: depl(MSIZE_TDOFS_VEC)
    aster_logical :: l_axis
!
! --------------------------------------------------------------------------------------------------
!
! --- Get HHO informations
!
    call elrefe_info(fami=fami, npg=npg)
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)
!
! --- Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym)
    faces_dofs = total_dofs-cbs
    gbs_cmp = gbs/(hhoCell%ndim*hhoCell%ndim)
    cbs_cmp = cbs/hhoCell%ndim
    nbsig = nbsigm()
!
! --- Type of finite element
!
    call hhoBasisCell%initialize(hhoCell)
    l_axis = lteatt('AXIS', 'OUI')
!
! --- Compute Operators
!
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, ASTER_FALSE, gradsym)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, ASTER_FALSE, gradsym)
    end if
!
    call readVector("PDEPLAR", total_dofs, depl)
!
! - Current time
    time = r8vide()
!
! ----- Material parameters
    call jevech('PMATERC', 'L', jvMater)
!
    call rcangm(hhoCell%ndim, hhoCell%barycenter, anglNaut)
!
    sigmVarc = 0.d0
!
! - Compute stresses from external state variables
    call sigtmc(fami, nbsig, npg, hhoCell%ndim, &
                time, zi(jvMater), anglNaut, &
                VARC_STRAIN_ALL, sigmVarc)
!
    sief = 0.d0
!
! ----- compute E = gradrec_sym * depl
    call gradsym%dot(depl, E_coeff)
!
! ------ Compute total stress
!
    do ipg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
        weight = hhoQuadCellRigi%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
        call hhoBasisCell%BSEval(coorpg(1:3), 0, &
                                 max(hhoData%grad_degree(), hhoData%cell_degree()), &
                                 BSCEval)
!
! ------ Eval deformations
        Eps = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_coeff)
!
        if (l_axis) then
            call hhoAddAxisGradSym(hhoCell, BSCEval, depl(faces_dofs+1:), &
                                   coorpg, cbs_cmp, Eps)
        end if
!
! ----- Compute elasticity matrix
        call dmatmc(fami, zi(jvMater), time, '+', &
                    ipg, 1, anglNaut, nbsig, dmat(1:nbsig, 1:nbsig))
!
! ----- Compute SIGM_ELGA
        sigma = 0.d0
        do i = 1, nbsig
            do j = 1, nbsig
                sigma(i) = sigma(i)+Eps(j)*dmat(i, j)
            end do
        end do
!
! ----- Compute SIEF_ELGA
!
        sief((ipg-1)*nbsig+1:ipg*nbsig) = sigma(1:nbsig)+sigmVarc((ipg-1)*nbsig+1:ipg*nbsig)
    end do
!
! ----- Set output vector
    call writeVector("PCONTRR", npg*nbsig, sief)
!
    call gradsym%free()
!
end subroutine
