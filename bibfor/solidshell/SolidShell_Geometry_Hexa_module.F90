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
! ==================================================================================================
!
! Module for geometry of HEXA cells for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Geometry_Hexa_module
! ==================================================================================================
    use SolidShell_type
    use SolidShell_Geometry_module
    use SolidShell_Debug_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: decoJacoMatrHexa, initGeomCellHexa, compAhmadFrame
    private :: compTMatrHexa, decoTMatrHexa, compTMatrElemHexa
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/SolidShell_type.h"
#include "blas/ddot.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! initGeomCellHexa
!
! Prepare geometric quantities for HEXA cell
!
! In  cellGeom         : general geometric properties of cell
! Out geomHexa         : geometric properties for HEXA cell
! In  disp             : diplacement to add at initial geometry
!
! --------------------------------------------------------------------------------------------------
    subroutine initGeomCellHexa(cellGeom, geomHexa, disp_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_GEOM_HEXA), intent(out) :: geomHexa
        real(kind=8), optional, intent(in) :: disp_(SSH_NBDOF_HEXA)
!   ------------------------------------------------------------------------------------------------
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('> initGeomCellHexa')
!
! - Init
        geomHexa%cellGeom = cellGeom
!
! - Set configuration
        geomHexa%geomCurr = geomHexa%cellGeom%geomInit
        if (present(disp_)) then
            geomHexa%geomCurr(1:SSH_NBDOFG_HEXA) = geomHexa%geomCurr(1:SSH_NBDOFG_HEXA)+disp_(1:&
                                                   &SSH_NBDOFG_HEXA)
        end if
!
! - Compute T matrix (matrix relating the covariant and cartesian frames) for HEXA cell
        call compTMatrHexa(geomHexa)
!
! - Decomposition of T matrix
        call decoTMatrHexa(geomHexa)
!
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('< initGeomCellHexa')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compTMatrHexa
!
! Compute T matrix (matrix relating the covariant and cartesian frames) for HEXA cell
!
! IO  geomHexa         : geometric properties for HEXA cell
!
! --------------------------------------------------------------------------------------------------
    subroutine compTMatrHexa(geomHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(inout) :: geomHexa
! - Local
        real(kind=8) :: JacInv0(3, 3)
!   ------------------------------------------------------------------------------------------------
!
        geomHexa%T0 = 0.d0
!
! - Get parameters
        JacInv0 = geomHexa%cellGeom%JacInv0
!
! - Compute T0 matrix
        geomHexa%T0(1, 1) = JacInv0(1, 1)*JacInv0(1, 1)
        geomHexa%T0(1, 2) = JacInv0(2, 1)*JacInv0(2, 1)
        geomHexa%T0(1, 3) = JacInv0(3, 1)*JacInv0(3, 1)
        geomHexa%T0(1, 4) = JacInv0(1, 1)*JacInv0(2, 1)
        geomHexa%T0(1, 5) = JacInv0(1, 1)*JacInv0(3, 1)
        geomHexa%T0(1, 6) = JacInv0(2, 1)*JacInv0(3, 1)
        geomHexa%T0(2, 1) = JacInv0(1, 2)*JacInv0(1, 2)
        geomHexa%T0(2, 2) = JacInv0(2, 2)*JacInv0(2, 2)
        geomHexa%T0(2, 3) = JacInv0(3, 2)*JacInv0(3, 2)
        geomHexa%T0(2, 4) = JacInv0(1, 2)*JacInv0(2, 2)
        geomHexa%T0(2, 5) = JacInv0(1, 2)*JacInv0(3, 2)
        geomHexa%T0(2, 6) = JacInv0(2, 2)*JacInv0(3, 2)
        geomHexa%T0(3, 1) = JacInv0(1, 3)*JacInv0(1, 3)
        geomHexa%T0(3, 2) = JacInv0(2, 3)*JacInv0(2, 3)
        geomHexa%T0(3, 3) = JacInv0(3, 3)*JacInv0(3, 3)
        geomHexa%T0(3, 4) = JacInv0(1, 3)*JacInv0(2, 3)
        geomHexa%T0(3, 5) = JacInv0(1, 3)*JacInv0(3, 3)
        geomHexa%T0(3, 6) = JacInv0(2, 3)*JacInv0(3, 3)
        geomHexa%T0(4, 1) = 2.d0*JacInv0(1, 1)*JacInv0(1, 2)
        geomHexa%T0(4, 2) = 2.d0*JacInv0(2, 1)*JacInv0(2, 2)
        geomHexa%T0(4, 3) = 2.d0*JacInv0(3, 1)*JacInv0(3, 2)
        geomHexa%T0(4, 4) = JacInv0(1, 1)*JacInv0(2, 2)+JacInv0(2, 1)*JacInv0(1, 2)
        geomHexa%T0(4, 5) = JacInv0(1, 1)*JacInv0(3, 2)+JacInv0(3, 1)*JacInv0(1, 2)
        geomHexa%T0(4, 6) = JacInv0(2, 2)*JacInv0(3, 1)+JacInv0(2, 1)*JacInv0(3, 2)
        geomHexa%T0(6, 1) = 2.d0*JacInv0(1, 2)*JacInv0(1, 3)
        geomHexa%T0(6, 2) = 2.d0*JacInv0(2, 2)*JacInv0(2, 3)
        geomHexa%T0(6, 3) = 2.d0*JacInv0(3, 2)*JacInv0(3, 3)
        geomHexa%T0(6, 4) = JacInv0(1, 3)*JacInv0(2, 2)+JacInv0(1, 2)*JacInv0(2, 3)
        geomHexa%T0(6, 5) = JacInv0(1, 2)*JacInv0(3, 3)+JacInv0(3, 2)*JacInv0(1, 3)
        geomHexa%T0(6, 6) = JacInv0(2, 2)*JacInv0(3, 3)+JacInv0(3, 2)*JacInv0(2, 3)
        geomHexa%T0(5, 1) = 2.d0*JacInv0(1, 1)*JacInv0(1, 3)
        geomHexa%T0(5, 2) = 2.d0*JacInv0(2, 1)*JacInv0(2, 3)
        geomHexa%T0(5, 3) = 2.d0*JacInv0(3, 1)*JacInv0(3, 3)
        geomHexa%T0(5, 4) = JacInv0(1, 3)*JacInv0(2, 1)+JacInv0(1, 1)*JacInv0(2, 3)
        geomHexa%T0(5, 5) = JacInv0(1, 3)*JacInv0(3, 1)+JacInv0(1, 1)*JacInv0(3, 3)
        geomHexa%T0(5, 6) = JacInv0(2, 3)*JacInv0(3, 1)+JacInv0(2, 1)*JacInv0(3, 3)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! decoTMatrHexa
!
! Decomposition of T matrix for HEXA cell: TXI, TETA, TZETA
!
! IO  geomHexa         : geometric properties for HEXA cell
!
! --------------------------------------------------------------------------------------------------
    subroutine decoTMatrHexa(geomHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(inout) :: geomHexa
! - Local
        real(kind=8) :: JacInv0(3, 3)
        real(kind=8) :: A(3, 3), B(3, 3), C(3, 3)
        real(kind=8) :: geomInitX(SSH_NBNODEG_HEXA), geomInitY(SSH_NBNODEG_HEXA)
        real(kind=8) :: geomInitZ(SSH_NBNODEG_HEXA)
        real(kind=8) :: JXI(3, 3), JETA(3, 3), JZETA(3, 3)
        real(kind=8) :: JXIJ0(3, 3), JETAJ0(3, 3), JZETAJ0(3, 3)
!   ------------------------------------------------------------------------------------------------
!
        geomHexa%TXI = 0.d0
        geomHexa%TETA = 0.d0
        geomHexa%TZETA = 0.d0
!
! - Get parameters
        JacInv0 = geomHexa%cellGeom%JacInv0
        geomInitX = geomHexa%cellGeom%geomInitX
        geomInitY = geomHexa%cellGeom%geomInitY
        geomInitZ = geomHexa%cellGeom%geomInitZ
!
! - Compute
        JXI = 0.d0
        JXI(1, 1) = 0.d0
        JXI(1, 2) = sum(geomInitX*hexaVectH1)
        JXI(1, 3) = sum(geomInitX*hexaVectH3)
        JXI(2, 1) = 0.d0
        JXI(2, 2) = sum(geomInitY*hexaVectH1)
        JXI(2, 3) = sum(geomInitY*hexaVectH3)
        JXI(3, 1) = 0.d0
        JXI(3, 2) = sum(geomInitZ*hexaVectH1)
        JXI(3, 3) = sum(geomInitZ*hexaVectH3)
!
        JETA = 0.d0
        JETA(1, 1) = sum(geomInitX*hexaVectH1)
        JETA(1, 2) = 0.d0
        JETA(1, 3) = sum(geomInitX*hexaVectH2)
        JETA(2, 1) = sum(geomInitY*hexaVectH1)
        JETA(2, 2) = 0.d0
        JETA(2, 3) = sum(geomInitY*hexaVectH2)
        JETA(3, 1) = sum(geomInitZ*hexaVectH1)
        JETA(3, 2) = 0.d0
        JETA(3, 3) = sum(geomInitZ*hexaVectH2)
!
        JZETA = 0.d0
        JZETA(1, 1) = sum(geomInitX*hexaVectH3)
        JZETA(1, 2) = sum(geomInitX*hexaVectH2)
        JZETA(1, 3) = 0.d0
        JZETA(2, 1) = sum(geomInitY*hexaVectH3)
        JZETA(2, 2) = sum(geomInitY*hexaVectH2)
        JZETA(2, 3) = 0.d0
        JZETA(3, 1) = sum(geomInitZ*hexaVectH3)
        JZETA(3, 2) = sum(geomInitZ*hexaVectH2)
        JZETA(3, 3) = 0.d0
!
        JXIJ0 = 0.d0
        JETAJ0 = 0.d0
        JZETAJ0 = 0.d0
        JXIJ0 = matmul(JXI, JacInv0)
        JETAJ0 = matmul(JETA, JacInv0)
        JZETAJ0 = matmul(JZETA, JacInv0)
!
        A = 0.d0
        B = 0.d0
        C = 0.d0
        A = -matmul(JacInv0, JXIJ0)
        B = -matmul(JacInv0, JETAJ0)
        C = -matmul(JacInv0, JZETAJ0)
!
        call compTMatrElemHexa(JacInv0, A, geomHexa%TXI)
        call compTMatrElemHexa(JacInv0, B, geomHexa%TETA)
        call compTMatrElemHexa(JacInv0, C, geomHexa%TZETA)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compTMatrElemHexa
!
! Compute elementary T matrix for HEXA cell
!
! In  JacInv0          : inverse of jacobian matrix at center of element on initial configuration
! In  A                : multiplicative matrix
! Out matrTElem        : elementary T matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compTMatrElemHexa(JacInv0, A, matrTElem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in) :: JacInv0(3, 3), A(3, 3)
        real(kind=8), intent(out) :: matrTElem(6, 6)
!   ------------------------------------------------------------------------------------------------
!
        matrTElem(1, 1) = 2.d0*JacInv0(1, 1)*A(1, 1)
        matrTElem(1, 2) = 2.d0*JacInv0(2, 1)*A(2, 1)
        matrTElem(1, 3) = 2.d0*JacInv0(3, 1)*A(3, 1)
        matrTElem(1, 4) = JacInv0(1, 1)*A(2, 1)+JacInv0(2, 1)*A(1, 1)
        matrTElem(1, 5) = JacInv0(1, 1)*A(3, 1)+JacInv0(3, 1)*A(1, 1)
        matrTElem(1, 6) = JacInv0(2, 1)*A(3, 1)+JacInv0(3, 1)*A(2, 1)
!
        matrTElem(2, 1) = 2.d0*JacInv0(1, 2)*A(1, 2)
        matrTElem(2, 2) = 2.d0*JacInv0(2, 2)*A(2, 2)
        matrTElem(2, 3) = 2.d0*JacInv0(3, 2)*A(3, 2)
        matrTElem(2, 4) = JacInv0(1, 2)*A(2, 2)+JacInv0(2, 2)*A(1, 2)
        matrTElem(2, 5) = JacInv0(1, 2)*A(3, 2)+JacInv0(3, 2)*A(1, 2)
        matrTElem(2, 6) = JacInv0(2, 2)*A(3, 2)+JacInv0(3, 2)*A(2, 2)
!
        matrTElem(3, 1) = 2.d0*JacInv0(1, 3)*A(1, 3)
        matrTElem(3, 2) = 2.d0*JacInv0(2, 3)*A(2, 3)
        matrTElem(3, 3) = 2.d0*JacInv0(3, 3)*A(3, 3)
        matrTElem(3, 4) = JacInv0(1, 3)*A(2, 3)+JacInv0(2, 3)*A(1, 3)
        matrTElem(3, 5) = JacInv0(1, 3)*A(3, 3)+JacInv0(3, 3)*A(1, 3)
        matrTElem(3, 6) = JacInv0(2, 3)*A(3, 3)+JacInv0(3, 3)*A(2, 3)
!
        matrTElem(4, 1) = 2.d0*(JacInv0(1, 1)*A(1, 2)+JacInv0(1, 2)*A(1, 1))
        matrTElem(4, 2) = 2.d0*(JacInv0(2, 1)*A(2, 2)+JacInv0(2, 2)*A(2, 1))
        matrTElem(4, 3) = 2.d0*(JacInv0(3, 1)*A(3, 2)+JacInv0(3, 2)*A(3, 1))
        matrTElem(4, 4) = ( &
                          JacInv0(1, 2)*A(2, 1)+JacInv0(2, 1)*A(1, 2))+(JacInv0(1, 1)*A(2, 2)+Ja&
                          &cInv0(2, 2)*A(1, 1) &
                          )
        matrTElem(4, 5) = ( &
                          JacInv0(1, 2)*A(3, 1)+JacInv0(3, 1)*A(1, 2))+(JacInv0(1, 1)*A(3, 2)+Ja&
                          &cInv0(3, 2)*A(1, 1) &
                          )
        matrTElem(4, 6) = ( &
                          JacInv0(2, 2)*A(3, 1)+JacInv0(3, 1)*A(2, 2))+(JacInv0(2, 1)*A(3, 2)+Ja&
                          &cInv0(3, 2)*A(2, 1) &
                          )
!
        matrTElem(5, 1) = 2.d0*(JacInv0(1, 1)*A(1, 3)+JacInv0(1, 3)*A(1, 1))
        matrTElem(5, 2) = 2.d0*(JacInv0(2, 1)*A(2, 3)+JacInv0(2, 3)*A(2, 1))
        matrTElem(5, 3) = 2.d0*(JacInv0(3, 1)*A(3, 3)+JacInv0(3, 3)*A(3, 1))
        matrTElem(5, 4) = ( &
                          JacInv0(1, 3)*A(2, 1)+JacInv0(2, 1)*A(1, 3))+(JacInv0(1, 1)*A(2, 3)+Ja&
                          &cInv0(2, 3)*A(1, 1) &
                          )
        matrTElem(5, 5) = ( &
                          JacInv0(1, 3)*A(3, 1)+JacInv0(3, 1)*A(1, 3))+(JacInv0(1, 1)*A(3, 3)+Ja&
                          &cInv0(3, 3)*A(1, 1) &
                          )
        matrTElem(5, 6) = ( &
                          JacInv0(2, 3)*A(3, 1)+JacInv0(3, 1)*A(2, 3))+(JacInv0(2, 1)*A(3, 3)+Ja&
                          &cInv0(3, 3)*A(2, 1) &
                          )
!
        matrTElem(6, 1) = 2.d0*(JacInv0(1, 2)*A(1, 3)+JacInv0(1, 3)*A(1, 2))
        matrTElem(6, 2) = 2.d0*(JacInv0(2, 2)*A(2, 3)+JacInv0(2, 3)*A(2, 2))
        matrTElem(6, 3) = 2.d0*(JacInv0(3, 2)*A(3, 3)+JacInv0(3, 3)*A(3, 2))
        matrTElem(6, 4) = ( &
                          JacInv0(1, 3)*A(2, 2)+JacInv0(2, 2)*A(1, 3))+(JacInv0(1, 2)*A(2, 3)+Ja&
                          &cInv0(2, 3)*A(1, 2) &
                          )
        matrTElem(6, 5) = ( &
                          JacInv0(1, 3)*A(3, 2)+JacInv0(3, 2)*A(1, 3))+(JacInv0(1, 2)*A(3, 3)+Ja&
                          &cInv0(3, 3)*A(1, 2) &
                          )
        matrTElem(6, 6) = ( &
                          JacInv0(2, 3)*A(3, 2)+JacInv0(3, 2)*A(2, 3))+(JacInv0(2, 2)*A(3, 3)+Ja&
                          &cInv0(3, 3)*A(2, 2) &
                          )
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! decoJacoMatrHexa
!
! Compute the decomposed derivatives of the jacobian in the covariant base for HEXA cell
!
! Decomposition of jacobian:
!   J  = J0 + xi.JXI + eta.JETA + zeta.JZETA +
!        xi.eta.JXIETA + eta.zeta.JETAZETA + xi.zeta.JXIZETA
!
! In  geomCurr         : current configuration of element
! Out J10              : decomposed derivatives of the jacobian in the covariant base
! Out J20              : decomposed derivatives of the jacobian in the covariant base
! Out J30              : decomposed derivatives of the jacobian in the covariant base
! Out J1ETA            : decomposed derivatives of the jacobian in the covariant base
! Out J1ZETA           : decomposed derivatives of the jacobian in the covariant base
! Out J1ETAZETA        : decomposed derivatives of the jacobian in the covariant base
! Out J2XI             : decomposed derivatives of the jacobian in the covariant base
! Out J2ZETA           : decomposed derivatives of the jacobian in the covariant base
! Out J2XIZETA         : decomposed derivatives of the jacobian in the covariant base
! Out J3XI             : decomposed derivatives of the jacobian in the covariant base
! Out J3ETA            : decomposed derivatives of the jacobian in the covariant base
! Out J3XIETA          : decomposed derivatives of the jacobian in the covariant base
!
! --------------------------------------------------------------------------------------------------
    subroutine decoJacoMatrHexa(geomCurr, J10, J1ETA, J1ZETA, J1ETAZETA, &
                                J20, J2XI, J2ZETA, J2XIZETA, J30, &
                                J3ETA, J3XI, J3XIETA)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in) :: geomCurr(SSH_NBDOFG_HEXA)
        real(kind=8), intent(out) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
        real(kind=8), intent(out) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
        real(kind=8), intent(out) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
! - Local
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom
        real(kind=8) :: XCurr(SSH_NBNODEG_HEXA), YCurr(SSH_NBNODEG_HEXA), ZCurr(SSH_NBNODEG_HEXA)
!   ------------------------------------------------------------------------------------------------
!
        J10 = 0.d0
        J1ETA = 0.d0
        J1ZETA = 0.d0
        J1ETAZETA = 0.d0
        J20 = 0.d0
        J2XI = 0.d0
        J2ZETA = 0.d0
        J2XIZETA = 0.d0
        J30 = 0.d0
        J3ETA = 0.d0
        J3XI = 0.d0
        J3XIETA = 0.d0
!
! - Get parameters
        do iNodeGeom = 1, nbNodeGeom
            XCurr(iNodeGeom) = geomCurr(3*(iNodeGeom-1)+1)
            YCurr(iNodeGeom) = geomCurr(3*(iNodeGeom-1)+2)
            ZCurr(iNodeGeom) = geomCurr(3*(iNodeGeom-1)+3)
        end do
!
! - Compute
        J10(1) = sum(hexaVectG1*XCurr)
        J10(2) = sum(hexaVectG1*YCurr)
        J10(3) = sum(hexaVectG1*ZCurr)
        J1ETA(1) = sum(hexaVectH1*XCurr)
        J1ETA(2) = sum(hexaVectH1*YCurr)
        J1ETA(3) = sum(hexaVectH1*ZCurr)
        J1ZETA(1) = sum(hexaVectH3*XCurr)
        J1ZETA(2) = sum(hexaVectH3*YCurr)
        J1ZETA(3) = sum(hexaVectH3*ZCurr)
        J1ETAZETA(1) = sum(hexaVectH4*XCurr)
        J1ETAZETA(2) = sum(hexaVectH4*YCurr)
        J1ETAZETA(3) = sum(hexaVectH4*ZCurr)
        J20(1) = sum(hexaVectG2*XCurr)
        J20(2) = sum(hexaVectG2*YCurr)
        J20(3) = sum(hexaVectG2*ZCurr)
        J2XI(1) = sum(hexaVectH1*XCurr)
        J2XI(2) = sum(hexaVectH1*YCurr)
        J2XI(3) = sum(hexaVectH1*ZCurr)
        J2ZETA(1) = sum(hexaVectH2*XCurr)
        J2ZETA(2) = sum(hexaVectH2*YCurr)
        J2ZETA(3) = sum(hexaVectH2*ZCurr)
        J2XIZETA(1) = sum(hexaVectH4*XCurr)
        J2XIZETA(2) = sum(hexaVectH4*YCurr)
        J2XIZETA(3) = sum(hexaVectH4*ZCurr)
        J30(1) = sum(hexaVectG3*XCurr)
        J30(2) = sum(hexaVectG3*YCurr)
        J30(3) = sum(hexaVectG3*ZCurr)
        J3ETA(1) = sum(hexaVectH2*XCurr)
        J3ETA(2) = sum(hexaVectH2*YCurr)
        J3ETA(3) = sum(hexaVectH2*ZCurr)
        J3XI(1) = sum(hexaVectH3*XCurr)
        J3XI(2) = sum(hexaVectH3*YCurr)
        J3XI(3) = sum(hexaVectH3*ZCurr)
        J3XIETA(1) = sum(hexaVectH4*XCurr)
        J3XIETA(2) = sum(hexaVectH4*YCurr)
        J3XIETA(3) = sum(hexaVectH4*ZCurr)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compAhmadFrame
!
! Compute quantities in Ahmad frame for pinch quantities
!
! In  cellGeom         : general geometric properties of cell
! Out area             : area
!
! --------------------------------------------------------------------------------------------------
    subroutine compAhmadFrame(cellGeom, area)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(out) :: area
! - Local
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom
        real(kind=8) :: XCurr(SSH_NBNODEG_HEXA), YCurr(SSH_NBNODEG_HEXA), ZCurr(SSH_NBNODEG_HEXA)
        real(kind=8) :: d, nX, nY, nZ, ps
        real(kind=8) :: U(3), V(3), W(3), vectAhmad2NonNorm(3)
        real(kind=8) :: vectAhmad1(3), vectAhmad2(3), vectAhmad3(3)
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
        area = 0.d0
!
! - Get parameters
        do iNodeGeom = 1, nbNodeGeom
            XCurr(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+1)
            YCurr(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+2)
            ZCurr(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+3)
        end do
!
! - Construct first axis
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        U(1) = ddot(b_n, XCurr, b_incx, hexaVectG1*8.d0, b_incy)
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        U(2) = ddot(b_n, YCurr, b_incx, hexaVectG1*8.d0, b_incy)
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        U(3) = ddot(b_n, ZCurr, b_incx, hexaVectG1*8.d0, b_incy)
        d = sqrt(U(1)*U(1)+U(2)*U(2)+U(3)*U(3))
        vectAhmad1(1) = U(1)/d
        vectAhmad1(2) = U(2)/d
        vectAhmad1(3) = U(3)/d
!
! - Construct second axis
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        V(1) = ddot(b_n, XCurr, b_incx, hexaVectG2*8.d0, b_incy)
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        V(2) = ddot(b_n, YCurr, b_incx, hexaVectG2*8.d0, b_incy)
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        V(3) = ddot(b_n, ZCurr, b_incx, hexaVectG2*8.d0, b_incy)
        vectAhmad2NonNorm = V
        ps = vectAhmad1(1)*V(1)+vectAhmad1(2)*V(2)+vectAhmad1(3)*V(3)
        V(1) = V(1)-ps*vectAhmad1(1)
        V(2) = V(2)-ps*vectAhmad1(2)
        V(3) = V(3)-ps*vectAhmad1(3)
        d = sqrt(V(1)*V(1)+V(2)*V(2)+V(3)*V(3))
        vectAhmad2(1) = V(1)/d
        vectAhmad2(2) = V(2)/d
        vectAhmad2(3) = V(3)/d
!
! - Construct third axis
        W(1) = vectAhmad1(2)*vectAhmad2(3)-vectAhmad1(3)*vectAhmad2(2)
        W(2) = vectAhmad1(3)*vectAhmad2(1)-vectAhmad1(1)*vectAhmad2(3)
        W(3) = vectAhmad1(1)*vectAhmad2(2)-vectAhmad1(2)*vectAhmad2(1)
        d = sqrt(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
        vectAhmad3(1) = W(1)/d
        vectAhmad3(2) = W(2)/d
        vectAhmad3(3) = W(3)/d
!
! - Area
        nX = U(2)*vectAhmad2NonNorm(3)-U(3)*vectAhmad2NonNorm(2)
        nY = U(3)*vectAhmad2NonNorm(1)-U(1)*vectAhmad2NonNorm(3)
        nZ = U(1)*vectAhmad2NonNorm(2)-U(2)*vectAhmad2NonNorm(1)
        area = sqrt(nX*nX+nY*nY+nZ*nZ)/16.0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Geometry_Hexa_module
