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
module coupling_quadrature_module
!
    use FE_quadrature_module
    use FE_topo_module
    use HHO_quadrature_module
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_penalisation_module.h"
#include "asterfort/jevech.h"
#include "asterfort/latrco.h"
#include "asterfort/lcptga.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmmjac.h"
#include "asterfort/mmnewd.h"
#include "asterfort/mmnonf.h"
#include "asterfort/projInsideCell.h"
#include "asterfort/reerel.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - Quadrature
!
! --------------------------------------------------------------------------------------------------
!
!
#define MAX_NB_INTE  8
#define PROJ_TOLE  1.0d-8

!
    private :: cplGetQuad, getQuadOrderFEM, cplProjFEM2HHO
    public :: cplGetQuadFEMHHO
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    integer(kind=8) function getQuadOrderFEM(typema)
!
        implicit none
!
        character(len=8), intent(in) :: typema
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        select case (typema)
        case ("SE2", "TR3", "QU4")
            getQuadOrderFEM = 1
        case ("SE3", "TR6", "TR7", "QU8", "QU9")
            getQuadOrderFEM = 2
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplGetQuad(feFace, max_order, FEQuad)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        integer(kind=8), intent(in) :: max_order
        type(FE_Quadrature), intent(out) :: FEQuad
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: nbPoinInte, jpair, iPoinInte, ndim
        integer(kind=8) :: nbPtGauss, nb_qp, i
        integer(kind=8) :: iTria, iGauss, nbTria, nbGauss
        real(kind=8) :: triaCoorSlav(2, 3), jaco, coorac(3)
        real(kind=8) :: gausWeightSlav(12), gausCoorSlav(2, 12)
        real(kind=8) :: shape_func(9), shape_dfunc(2, 9)
        character(len=8) :: elga_fami
        real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
!
        call jevech('PPAIRR', 'L', jpair)
        nbPoinInte = int(zr(jpair-1+2))
        ASSERT(nbPoinInte .le. MAX_NB_INTE)
        poinInteSlav = 0.d0
        do iPoinInte = 1, nbPoinInte
            poinInteSlav(1, iPoinInte) = zr(jpair-1+2+iPoinInte)
            poinInteSlav(2, iPoinInte) = zr(jpair-1+10+iPoinInte)
        end do
!
        ndim = feFace%ndim+1
!
! - Triangulation of convex polygon defined by intersection points
        if (ndim .eq. 3) then
            if (nbPoinInte == 3) then
                nbTria = 1
            else
                nbTria = nbPoinInte
            end if
            if (max_order <= 5) then
                ! order 5 - by triangle
                elga_fami = 'FPG4'
                nbPtGauss = 4
            else if (max_order <= 7) then
                ! order 7 - by triangle
                elga_fami = 'FPG7'
                nbPtGauss = 7
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (ndim .eq. 2) then
            nbTria = 1
!
            if (max_order <= 5) then
                ! order 5
                elga_fami = 'FPG3'
                nbPtGauss = 3
            else if (max_order <= 7) then
                ! order 7
                elga_fami = 'FPG4'
                nbPtGauss = 4
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
        ASSERT(nbTria*nbPtGauss <= 64)
        FEQuad%nbQuadPoints = nbTria*nbPtGauss
        nb_qp = 0
!
! - Loop on triangles
        do iTria = 1, nbTria
! ----- Coordinates of current triangle (slave)
            triaCoorSlav = 0.d0
            if (ndim .eq. 3) then
                call latrco(iTria, nbPoinInte, poinInteSlav, triaCoorSlav)
            elseif (ndim .eq. 2) then
                ASSERT(nbPoinInte .eq. 2)
                triaCoorSlav(1:2, 1:2) = poinInteSlav(1:2, 1:2)
            end if

! ----- Get integration points for slave element
            call lcptga(ndim, triaCoorSlav, elga_fami, &
                        nbGauss, gausCoorSlav, gausWeightSlav)

! ----- Loop on integration points
            do iGauss = 1, nbGauss
                nb_qp = nb_qp+1
                FEQuad%points_param(1:2, nb_qp) = gausCoorSlav(1:2, iGauss)
                FEQuad%weights_param(nb_qp) = gausWeightSlav(iGauss)
!
! ------------- Get shape functions and first derivative only (for perf)
                call mmnonf(ndim, feFace%nbnodes, feFace%typemas, &
                            gausCoorSlav(1, iGauss), gausCoorSlav(2, iGauss), &
                            shape_func)
                call mmdonf(ndim, feFace%nbnodes, feFace%typemas, &
                            gausCoorSlav(1, iGauss), gausCoorSlav(2, iGauss), &
                            shape_dfunc)
!
                coorac = 0.d0
                do i = 1, feFace%nbnodes
                    coorac(1:3) = coorac(1:3)+feFace%coorno(1:3, i)*shape_func(i)
                end do
!
! ------------- Compute jacobian
                call mmmjac(ASTER_FALSE, feFace%nbnodes, ndim, &
                            feFace%typemas, feFace%coorno, &
                            shape_func, shape_dfunc, jaco)
!
                FEQuad%points(1:3, nb_qp) = coorac
                FEQuad%weights(nb_qp) = abs(jaco)*gausWeightSlav(iGauss)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplProjFEM2HHO(feFace, hhoFace, FEQuad, hhoQuad)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        type(HHO_Face), intent(in) :: hhoFace
        type(FE_Quadrature), intent(in) :: FEQuad
        type(HHO_Quadrature), intent(out) :: hhoQuad
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: ndim, ipg, iret
        real(kind=8) :: norm_slav(3), tau_slav(3, 2), coor_qp_sl(2)
        real(kind=8) :: coor_qp_sl_re(3), ksi_line(2)
        real(kind=8) :: tau1_mast(3), tau2_mast(3)
        character(len=8) :: type_mast
!
        hhoQuad%nbQuadPoints = FEQuad%nbQuadPoints
        ndim = feFace%ndim+1
!
        call CellNameL2S(hhoFace%typema, type_mast)
!
        do ipg = 1, hhoQuad%nbQuadPoints
            hhoQuad%weights(ipg) = FEQuad%weights(ipg)
            coor_qp_sl = FEQuad%points_param(1:2, ipg)
!
! ------ Compute outward slave normal (pairing configuration)
!
            call apnorm(feface%nbnodes, feface%typemas, ndim, feface%coorno, &
                        coor_qp_sl(1), coor_qp_sl(2), norm_slav, tau_slav(1:3, 1), tau_slav(1:3, 2))
!
! ----- Return in real slave space (pairing configuration)
!
            coor_qp_sl_re = FEQuad%points(1:3, ipg)
!
! ----- Projection on master element
            call mmnewd(type_mast, hhoFace%nbnodes, ndim, hhoFace%coorno, &
                        coor_qp_sl_re, 100, PROJ_TOLE, norm_slav, ksi_line(1), &
                        ksi_line(2), tau1_mast, tau2_mast, iret)
            ASSERT(iret == 0)
!
! ----- Check that projected node is inside cell
!
            call projInsideCell(1e-3, ndim, type_mast, ksi_line, iret)
            ASSERT(iret == 0)
!
            hhoQuad%points_param(1:2, ipg) = ksi_line
            call reerel(type_mast, hhoFace%nbnodes, 3, hhoFace%coorno, &
                        hhoQuad%points_param(1:3, ipg), &
                        hhoQuad%points(1:3, ipg))
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplGetQuadFEMHHO(feFace, hhoFace, hhoData, FEQuad, hhoQuad)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Data), intent(in) :: hhoData
        type(FE_Quadrature), intent(out) :: FEQuad
        type(HHO_Quadrature), intent(out) :: hhoQuad
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: max_order
!
        max_order = 2*(max(hhoData%face_degree(), getQuadOrderFEM(feFace%typemas))+1)
!
        call cplGetQuad(feFace, max_order, FEQuad)
        call cplProjFEM2HHO(feFace, hhoFace, FEQuad, hhoQuad)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
