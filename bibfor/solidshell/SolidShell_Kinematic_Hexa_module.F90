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
! Module for kinematic of HEXA cells for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Kinematic_Hexa_module
! ==================================================================================================
    use SolidShell_type
    use SolidShell_Geometry_module
    use SolidShell_Geometry_Hexa_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: compBCovaMatrHexa, compBCartMatrHexa, compBCartEASMatrHexa, &
              compBMatrHexa, compEpsiHexa, &
              compECovaMatrHexa, compEpsgHexa, &
              compEpslHexa
    public :: compGCovaMatrHexa, compGCovaSMatrHexa
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/SolidShell_type.h"
#include "asterc/r8miem.h"
#include "asterfort/diagp3.h"
#include "asterfort/tnsvec.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! compECovaMatrHexa
!
! Compute decomposition of Green-Lagrange strains in covariant basis for HEXA cell
!
! In  cellGeom         : general geometric properties of cell
! In  disp             : displacements of element
! IO  epsg             : Green-Lagrange strains
!
! --------------------------------------------------------------------------------------------------
    subroutine compECovaMatrHexa(cellGeom, disp, epsg)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_CELL_GEOM), intent(in)     :: cellGeom
        real(kind=8), intent(in)            :: disp(SSH_NBDOF_HEXA)
        type(SSH_EPSG_HEXA), intent(inout)  :: epsg
! - Local
        integer(kind=8) :: EH, AD, JM
        real(kind=8) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
        real(kind=8) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
        real(kind=8) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
        real(kind=8) :: J1(3), J2(3), J3(3)
        real(kind=8) :: XI(3)
        real(kind=8) :: D10(3), D1ETA(3), D1ZETA(3), D1ETAZETA(3)
        real(kind=8) :: D20(3), D2XI(3), D2ZETA(3), D2XIZETA(3)
        real(kind=8) :: D30(3), D3ETA(3), D3XI(3), D3XIETA(3)
        real(kind=8) :: D1(3), D2(3), D3(3)
!   ------------------------------------------------------------------------------------------------
!
        epsg%ECova0 = 0.d0
        epsg%ECovaZETA = 0.d0
        epsg%ECovaZETAZETA = 0.d0
        epsg%ECovaXI = 0.d0
        epsg%ECovaETA = 0.d0
        epsg%ECovaXIZETA = 0.d0
        epsg%ECovaETAZETA = 0.d0

! - Compute the decomposed derivatives of the jacobian in the covariant base
        call decoJacoMatrHexa(cellGeom%geomInit, &
                              J10, J1ETA, J1ZETA, J1ETAZETA, &
                              J20, J2XI, J2ZETA, J2XIZETA, &
                              J30, J3ETA, J3XI, J3XIETA)

! - Compute the decomposed derivatives of the displacements in the covariant base
        call decoJacoMatrHexa(disp(1:SSH_NBDOFG_HEXA), &
                              D10, D1ETA, D1ZETA, D1ETAZETA, &
                              D20, D2XI, D2ZETA, D2XIZETA, &
                              D30, D3ETA, D3XI, D3XIETA)

! - Compute ECova0
        epsg%ECova0(1) = sum(J10*D10)+0.5d0*sum(D10*D10)
        epsg%ECova0(2) = sum(J20*D20)+0.5d0*sum(D20*D20)
        epsg%ECova0(4) = sum(J10*D20)+sum(J20*D10)+sum(D10*D20)

! - Compute ECovaZETA
        epsg%ECovaZETA(1) = sum(J10*D1ZETA)+sum(D10*J1ZETA)+sum(D10*D1ZETA)
        epsg%ECovaZETA(2) = sum(J20*D2ZETA)+sum(D20*J2ZETA)+sum(D20*D2ZETA)
        epsg%ECovaZETA(4) = sum(J10*D2ZETA)+sum(D20*J1ZETA)+sum(J20*D1ZETA)+ &
                            sum(D10*J2ZETA)+sum(D10*D2ZETA)+sum(D20*D1ZETA)

! - Compute ECovaZETAZETA
        epsg%ECovaZETAZETA(1) = sum(J1ZETA*D1ZETA)+0.5d0*sum(D1ZETA*D1ZETA)
        epsg%ECovaZETAZETA(2) = sum(J2ZETA*D2ZETA)+0.5d0*sum(D2ZETA*D2ZETA)
        epsg%ECovaZETAZETA(4) = sum(J1ZETA*D2ZETA)+sum(J2ZETA*D1ZETA)+ &
                                sum(D2ZETA*D1ZETA)

! - Compute ECovaXI (only for stabilization)
        epsg%ECovaXI(2) = sum(J20*D2XI)+sum(J2XI*D20)+sum(D2XI*D20)
        epsg%ECovaXI(4) = sum(J10*D2XI)+sum(J2XI*D10)+sum(D2XI*D10)

! - Compute ECovaETA (only for stabilization)
        epsg%ECovaETA(1) = sum(J10*D1ETA)+sum(D10*J1ETA)+sum(D10*D1ETA)
        epsg%ECovaETA(4) = sum(J20*D1ETA)+sum(D20*J1ETA)+sum(D20*D1ETA)

! - Compute ECovaXIZETA (only for stabilization)
        epsg%ECovaXIZETA(2) = sum(J20*D2XIZETA)+sum(J2XI*D2ZETA)+sum(J2ZETA*D2XI)+ &
                              sum(D20*J2XIZETA)+sum(D20*D2XIZETA)+sum(D2XI*D2ZETA)
        epsg%ECovaXIZETA(4) = sum(J10*D2XIZETA)+sum(J2XI*D1ZETA)+sum(J1ZETA*D2XI)+ &
                              sum(D10*J2XIZETA)+sum(D10*D2XIZETA)+sum(D2XI*D1ZETA)

! - Compute ECovaETAZETA (only for stabilization)
        epsg%ECovaETAZETA(1) = sum(J10*D1ETAZETA)+sum(J1ETA*D1ZETA)+sum(J1ZETA*D1ETA)+ &
                               sum(D10*J1ETAZETA)+sum(D10*D1ETAZETA)+sum(D1ETA*D1ZETA)
        epsg%ECovaETAZETA(4) = sum(J20*D1ETAZETA)+sum(J1ETA*D2ZETA)+sum(J2ZETA*D1ETA)+ &
                               sum(D20*J1ETAZETA)+sum(D20*D1ETAZETA)+sum(D1ETA*D2ZETA)

! - Assumed strains
        do AD = 1, 4
            XI = hexaQuadXIAD(AD, 1:3)
            J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
            D3 = SSH_JACO_J3(D30, D3ETA, D3XI, D3XIETA, XI(1), XI(2))
            epsg%ECova0(3) = epsg%ECova0(3)+ &
                             (sum(J3*D3)+0.5d0*sum(D3*D3))/4.d0
            epsg%ECovaXI(3) = epsg%ECovaXI(3)+ &
                              XI(1)*(sum(J3*D3)+0.5d0*sum(D3*D3))/4.d0
            epsg%ECovaETA(3) = epsg%ECovaETA(3)+ &
                               XI(2)*(sum(J3*D3)+0.5d0*sum(D3*D3))/4.d0
        end do
        do EH = 1, 4
            XI = hexaQuadXIEH(EH, 1:3)
            J2 = SSH_JACO_J2(J20, J2XI, J2ZETA, J2XIZETA, XI(1), XI(3))
            J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
            D2 = SSH_JACO_J2(D20, D2XI, D2ZETA, D2XIZETA, XI(1), XI(3))
            D3 = SSH_JACO_J3(D30, D3ETA, D3XI, D3XIETA, XI(1), XI(2))
            epsg%ECova0(6) = epsg%ECova0(6)+ &
                             (sum(J3*D2)+sum(J2*D3)+sum(D3*D2))/4.d0
            epsg%ECovaZETA(6) = epsg%ECovaZETA(6)+ &
                                XI(3)*(sum(J3*D2)+sum(J2*D3)+sum(D3*D2))/4.d0
            epsg%ECovaXI(6) = epsg%ECovaXI(6)+ &
                              XI(1)*(sum(J2*D3)+sum(J3*D2)+sum(D2*D3))/4.d0
            epsg%ECovaXIZETA(6) = epsg%ECovaXIZETA(6)+ &
                                  XI(1)*XI(3)*(sum(J2*D3)+sum(J3*D2)+sum(D2*D3))/4.d0
        end do
        do JM = 1, 4
            XI = hexaQuadXIJM(JM, 1:3)
            J1 = SSH_JACO_J1(J10, J1ETA, J1ZETA, J1ETAZETA, XI(2), XI(3))
            J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
            D1 = SSH_JACO_J1(D10, D1ETA, D1ZETA, D1ETAZETA, XI(2), XI(3))
            D3 = SSH_JACO_J3(D30, D3ETA, D3XI, D3XIETA, XI(1), XI(2))
            epsg%ECova0(5) = epsg%ECova0(5)+ &
                             (sum(J3*D1)+sum(J1*D3)+sum(D3*D1))/4.d0
            epsg%ECovaZETA(5) = epsg%ECovaZETA(5)+ &
                                XI(3)*(sum(J3*D1)+sum(J1*D3)+sum(D3*D1))/4.d0
            epsg%ECovaETA(5) = epsg%ECovaETA(5)+ &
                               XI(2)*(sum(J1*D3)+sum(J3*D1)+sum(D1*D3))/4.d0
            epsg%ECovaETAZETA(5) = epsg%ECovaETAZETA(5)+ &
                                   XI(2)*XI(3)*(sum(J1*D3)+sum(J3*D1)+sum(D1*D3))/4.d0
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compBCovaMatrHexa
!
! Compute gradient matrix in covariant basis for HEXA cell
!
! In  geomHexa         : geometric properties for HEXA cell
! IO  kineHexa         : kinematic quantities (HEXA cell)
!
! --------------------------------------------------------------------------------------------------
    subroutine compBCovaMatrHexa(geomHexa, kineHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(in) :: geomHexa
        type(SSH_KINE_HEXA), intent(inout) :: kineHexa
! - Local
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom, EH, AD, JM, C1, C2
        real(kind=8) :: aux13(3), aux23(3), aux12(3), aux11(3), aux22(3), aux33(3)
        real(kind=8) :: dN_dXsi(SSH_NBNODEG_HEXA), dN_dEta(SSH_NBNODEG_HEXA)
        real(kind=8) :: dN_dZeta(SSH_NBNODEG_HEXA)
        real(kind=8) :: J10(3), J1ETA(3), J1ZETA(3), J1ETAZETA(3)
        real(kind=8) :: J20(3), J2XI(3), J2ZETA(3), J2XIZETA(3)
        real(kind=8) :: J30(3), J3ETA(3), J3XI(3), J3XIETA(3)
        real(kind=8) :: J1(3), J2(3), J3(3)
        real(kind=8) :: XI(3)
!   ------------------------------------------------------------------------------------------------
!
        kineHexa%BCova0 = 0.d0
        kineHexa%BCovaZETA = 0.d0
        kineHexa%BCovaZETAZETA = 0.d0
        kineHexa%BCovaXI = 0.d0
        kineHexa%BCovaETA = 0.d0
        kineHexa%BCovaETAZETA = 0.d0
        kineHexa%BCovaXIZETA = 0.d0

! - Compute the decomposed derivatives of the jacobian in the covariant base
        call decoJacoMatrHexa(geomHexa%geomCurr, &
                              J10, J1ETA, J1ZETA, J1ETAZETA, &
                              J20, J2XI, J2ZETA, J2XIZETA, &
                              J30, J3ETA, J3XI, J3XIETA)

! - Compute BCova0
        do iNodeGeom = 1, nbNodeGeom
            aux33 = 0.d0
            aux23 = 0.d0
            aux13 = 0.d0
            do AD = 1, 4
                XI = hexaQuadXIAD(AD, 1:3)
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux33 = aux33+dN_dZeta(iNodeGeom)*J3/4.d0
            end do
            do EH = 1, 4
                XI = hexaQuadXIEH(EH, 1:3)
                J2 = SSH_JACO_J2(J20, J2XI, J2ZETA, J2XIZETA, XI(1), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+(dN_dEta(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J2)/4.d0
            end do
            do JM = 1, 4
                XI = hexaQuadXIJM(JM, 1:3)
                J1 = SSH_JACO_J1(J10, J1ETA, J1ZETA, J1ETAZETA, XI(2), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+(dN_dXsi(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J1)/4.d0
            end do
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            kineHexa%BCova0(1, C1:C2) = hexaVectG1(iNodeGeom)*J10
            kineHexa%BCova0(2, C1:C2) = hexaVectG2(iNodeGeom)*J20
            kineHexa%BCova0(3, C1:C2) = aux33
            kineHexa%BCova0(4, C1:C2) = hexaVectG1(iNodeGeom)*J20+hexaVectG2(iNodeGeom)*J10
            kineHexa%BCova0(5, C1:C2) = aux13
            kineHexa%BCova0(6, C1:C2) = aux23
        end do

! - Compute BCovaZETA
        do iNodeGeom = 1, nbNodeGeom
            aux23 = 0.d0
            aux13 = 0.d0
            do EH = 1, 4
                XI = hexaQuadXIEH(EH, 1:3)
                J2 = SSH_JACO_J2(J20, J2XI, J2ZETA, J2XIZETA, XI(1), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+XI(3)*(dN_dEta(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J2)/4.d0
            end do
            do JM = 1, 4
                XI = hexaQuadXIJM(JM, 1:3)
                J1 = SSH_JACO_J1(J10, J1ETA, J1ZETA, J1ETAZETA, XI(2), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+XI(3)*(dN_dXsi(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J1)/4.d0
            end do
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            aux12 = hexaVectH2(iNodeGeom)*J10+ &
                    hexaVectH3(iNodeGeom)*J20+ &
                    hexaVectG2(iNodeGeom)*J1ZETA+ &
                    hexaVectG1(iNodeGeom)*J2ZETA
            kineHexa%BCovaZETA(1, C1:C2) = hexaVectH3(iNodeGeom)*J10+hexaVectG1(iNodeGeom)*J1ZETA
            kineHexa%BCovaZETA(2, C1:C2) = hexaVectH2(iNodeGeom)*J20+hexaVectG2(iNodeGeom)*J2ZETA
            kineHexa%BCovaZETA(3, C1:C2) = 0.d0
            kineHexa%BCovaZETA(4, C1:C2) = aux12
            kineHexa%BCovaZETA(5, C1:C2) = aux13
            kineHexa%BCovaZETA(6, C1:C2) = aux23
        end do

! - Compute BCovaZETAZETA
        do iNodeGeom = 1, nbNodeGeom
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            kineHexa%BCovaZETAZETA(1, C1:C2) = hexaVectH3(iNodeGeom)*J1ZETA
            kineHexa%BCovaZETAZETA(2, C1:C2) = hexaVectH2(iNodeGeom)*J2ZETA
            kineHexa%BCovaZETAZETA(3, C1:C2) = 0.d0
            kineHexa%BCovaZETAZETA(4, C1:C2) = hexaVectH2(iNodeGeom)* &
                                               J1ZETA+hexaVectH3(iNodeGeom)*J2ZETA
            kineHexa%BCovaZETAZETA(5, C1:C2) = 0.d0
            kineHexa%BCovaZETAZETA(6, C1:C2) = 0.d0
        end do

! - Compute BCovaXI
        do iNodeGeom = 1, nbNodeGeom
            aux33 = 0.d0
            aux23 = 0.d0
            do AD = 1, 4
                XI = hexaQuadXIAD(AD, 1:3)
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux33 = aux33+XI(1)*(dN_dZeta(iNodeGeom)*J3)/4.d0
            end do
            do EH = 1, 4
                XI = hexaQuadXIEH(EH, 1:3)
                J2 = SSH_JACO_J2(J20, J2XI, J2ZETA, J2XIZETA, XI(1), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+XI(1)*(dN_dEta(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J2)/4.d0
            end do
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            kineHexa%BCovaXI(1, C1:C2) = 0.d0
            kineHexa%BCovaXI(2, C1:C2) = hexaVectH1(iNodeGeom)*J20+hexaVectG2(iNodeGeom)*J2XI
            kineHexa%BCovaXI(3, C1:C2) = aux33
            kineHexa%BCovaXI(4, C1:C2) = hexaVectH1(iNodeGeom)*J10+hexaVectG1(iNodeGeom)*J2XI
            kineHexa%BCovaXI(5, C1:C2) = 0.d0
            kineHexa%BCovaXI(6, C1:C2) = aux23
        end do

! - Compute BCovaETA
        do iNodeGeom = 1, nbNodeGeom
            aux33 = 0.d0
            aux13 = 0.d0
            do AD = 1, 4
                XI = hexaQuadXIAD(AD, 1:3)
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux33 = aux33+XI(2)*(dN_dZeta(iNodeGeom)*J3)/4.d0
            end do
            do JM = 1, 4
                XI = hexaQuadXIJM(JM, 1:3)
                J1 = SSH_JACO_J1(J10, J1ETA, J1ZETA, J1ETAZETA, XI(2), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+XI(2)*(dN_dXsi(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J1)/4.d0
            end do
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            kineHexa%BCovaETA(1, C1:C2) = hexaVectH1(iNodeGeom)*J10+hexaVectG1(iNodeGeom)*J1ETA
            kineHexa%BCovaETA(2, C1:C2) = 0.d0
            kineHexa%BCovaETA(3, C1:C2) = aux33
            kineHexa%BCovaETA(4, C1:C2) = hexaVectH1(iNodeGeom)*J20+hexaVectG2(iNodeGeom)*J1ETA
            kineHexa%BCovaETA(5, C1:C2) = aux13
            kineHexa%BCovaETA(6, C1:C2) = 0.d0
        end do

! - Compute BCovaETAZETA
        do iNodeGeom = 1, nbNodeGeom
            aux11 = hexaVectH4(iNodeGeom)*J10+ &
                    hexaVectH3(iNodeGeom)*J1ETA+ &
                    hexaVectH1(iNodeGeom)*J1ZETA+ &
                    hexaVectG1(iNodeGeom)*J1ETAZETA
            aux12 = hexaVectH4(iNodeGeom)*J20+ &
                    hexaVectH2(iNodeGeom)*J1ETA+ &
                    hexaVectH1(iNodeGeom)*J2ZETA+ &
                    hexaVectG2(iNodeGeom)*J1ETAZETA
            aux13 = 0.d0
            do JM = 1, 4
                XI = hexaQuadXIJM(JM, 1:3)
                J1 = SSH_JACO_J1(J10, J1ETA, J1ZETA, J1ETAZETA, XI(2), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux13 = aux13+XI(2)*XI(3)*(dN_dXsi(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J1)/4.d0
            end do
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            kineHexa%BCovaETAZETA(1, C1:C2) = aux11
            kineHexa%BCovaETAZETA(2, C1:C2) = 0.d0
            kineHexa%BCovaETAZETA(3, C1:C2) = 0.d0
            kineHexa%BCovaETAZETA(4, C1:C2) = aux12
            kineHexa%BCovaETAZETA(5, C1:C2) = aux13
            kineHexa%BCovaETAZETA(6, C1:C2) = 0.d0
        end do

! - Compute BCovaXIZETA
        do iNodeGeom = 1, nbNodeGeom
            aux22 = hexaVectH4(iNodeGeom)*J20+ &
                    hexaVectH2(iNodeGeom)*J2XI+ &
                    hexaVectH1(iNodeGeom)*J2ZETA+ &
                    hexaVectG2(iNodeGeom)*J2XIZETA
            aux12 = hexaVectH4(iNodeGeom)*J10+ &
                    hexaVectH3(iNodeGeom)*J2XI+ &
                    hexaVectH1(iNodeGeom)*J1ZETA+ &
                    hexaVectG1(iNodeGeom)*J2XIZETA
            aux23 = 0.d0
            do EH = 1, 4
                XI = hexaQuadXIEH(EH, 1:3)
                J2 = SSH_JACO_J2(J20, J2XI, J2ZETA, J2XIZETA, XI(1), XI(3))
                J3 = SSH_JACO_J3(J30, J3ETA, J3XI, J3XIETA, XI(1), XI(2))
                call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
                aux23 = aux23+XI(1)*XI(3)*(dN_dEta(iNodeGeom)*J3+dN_dZeta(iNodeGeom)*J2)/4.d0
            end do
            C1 = 3*(iNodeGeom-1)+1
            C2 = 3*(iNodeGeom-1)+3
            kineHexa%BCovaXIZETA(1, C1:C2) = 0.d0
            kineHexa%BCovaXIZETA(2, C1:C2) = aux22
            kineHexa%BCovaXIZETA(3, C1:C2) = 0.d0
            kineHexa%BCovaXIZETA(4, C1:C2) = aux12
            kineHexa%BCovaXIZETA(5, C1:C2) = 0.d0
            kineHexa%BCovaXIZETA(6, C1:C2) = aux23
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compBCartMatrHexa
!
! Compute gradient matrix in cartesian frame for HEXA cell
!
! In  geomHexa         : geometric properties for HEXA cell
! IO  kineHexa         : kinematic quantities for HEXA cell
!
! --------------------------------------------------------------------------------------------------
    subroutine compBCartMatrHexa(geomHexa, kineHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(in)   :: geomHexa
        type(SSH_KINE_HEXA), intent(inout) :: kineHexa
!   ------------------------------------------------------------------------------------------------
!
        kineHexa%BCart0 = 0.d0
        kineHexa%BCartZ = 0.d0
        kineHexa%BCartZZ = 0.d0
        kineHexa%BCartX = 0.d0
        kineHexa%BCartY = 0.d0
        kineHexa%BCartYZ = 0.d0
        kineHexa%BCartXZ = 0.d0

        kineHexa%BCart0 = matmul(geomHexa%T0, kineHexa%BCova0)
        kineHexa%BCartZ = matmul(geomHexa%T0, kineHexa%BCovaZETA)+ &
                          matmul(geomHexa%TZETA, kineHexa%BCova0)
        kineHexa%BCartZZ = matmul(geomHexa%T0, kineHexa%BCovaZETAZETA)+ &
                           matmul(geomHexa%TZETA, kineHexa%BCovaZETA)
        kineHexa%BCartX = matmul(geomHexa%T0, kineHexa%BCovaXI)+ &
                          matmul(geomHexa%TXI, kineHexa%BCova0)
        kineHexa%BCartY = matmul(geomHexa%T0, kineHexa%BCovaETA)+ &
                          matmul(geomHexa%TETA, kineHexa%BCova0)
        kineHexa%BCartYZ = matmul(geomHexa%T0, kineHexa%BCovaETAZETA)+ &
                           matmul(geomHexa%TETA, kineHexa%BCovaZETA)+ &
                           matmul(geomHexa%TZETA, kineHexa%BCovaETA)
        kineHexa%BCartXZ = matmul(geomHexa%T0, kineHexa%BCovaXIZETA)+ &
                           matmul(geomHexa%TXI, kineHexa%BCovaZETA)+ &
                           matmul(geomHexa%TZETA, kineHexa%BCovaXI)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compBCartEASMatrHexa
!
! Compute EAS B matrix in cartesian frame for HEXA cell
!
! In  zeta             : out-of-plane parametric component
! In  geomHexa         : geometric properties for HEXA cell
! IO  kineHexa         : kinematic quantities (HEXA cell)
!
! --------------------------------------------------------------------------------------------------
    subroutine compBCartEASMatrHexa(zeta, geomHexa, kineHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in)           :: zeta
        real(kind=8)                       :: demih
        type(SSH_GEOM_HEXA), intent(in)    :: geomHexa
        type(SSH_KINE_HEXA), intent(inout) :: kineHexa
!   ------------------------------------------------------------------------------------------------
!
        kineHexa%BCartEAS = 0.d0
        demih = geomHexa%cellGeom%h0*0.5d0
        kineHexa%BCovaEAS = (/0.d0, 0.d0, -2.d0*zeta*demih, 0.d0, 0.d0, 0.d0/)
        kineHexa%BCartEAS = matmul(geomHexa%T0, kineHexa%BCovaEAS)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compBMatrHexa
!
! Compute B matrix for HEXA cell
!
! In  zeta             : out-of-plane parametric component
! IO  kineHexa         : kinematic quantities (HEXA cell)
!
! --------------------------------------------------------------------------------------------------
    subroutine compBMatrHexa(zeta, kineHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in) :: zeta
        type(SSH_KINE_HEXA), intent(inout) :: kineHexa
!   ------------------------------------------------------------------------------------------------
!
        kineHexa%B = 0.d0
        kineHexa%B(1:SSH_SIZE_TENS, 1:SSH_NBDOFG_HEXA) = kineHexa%BCart0+ &
                                                         zeta*kineHexa%BCartZ+ &
                                                         zeta*zeta*kineHexa%BCartZZ
        kineHexa%B(1:SSH_SIZE_TENS, SSH_NBDOF_HEXA) = kineHexa%BCartEAS
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsiHexa
!
! Compute small strains at Gauss points for HEXA cell
!
! In  kineHexa         : kinematic quantities (HEXA cell)
! In  disp             : displacements
! Out epsi             : small strains
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpsiHexa(kineHexa, disp, epsi)
! - Parameters
        type(SSH_KINE_HEXA), intent(in) :: kineHexa
        real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
        real(kind=8), intent(out)       :: epsi(SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
        epsi = 0.d0
        epsi = matmul(kineHexa%B(1:SSH_SIZE_TENS, 1:SSH_NBDOF_HEXA), disp(1:SSH_NBDOF_HEXA))
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsgHexa
!
! Compute Green-Lagrange strains in cartesian basis for HEXA cell
!
! In  zeta             : out-of-plane parametric component
! In  geomHexa         : geometric properties for HEXA cell
! IO  epsg             : Green-Lagrange strains
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpsgHexa(zeta, geomHexa, epsg)
! - Parameters
        real(kind=8), intent(in)           :: zeta
        type(SSH_GEOM_HEXA), intent(in)    :: geomHexa
        type(SSH_EPSG_HEXA), intent(inout) :: epsg
!   ------------------------------------------------------------------------------------------------
!
        epsg%vale = 0.d0
        epsg%vale = matmul(geomHexa%T0, epsg%ECova0)+ &
                    (matmul(geomHexa%T0, epsg%ECovaZETA)+ &
                     matmul(geomHexa%TZETA, epsg%ECova0))*zeta+ &
                    (matmul(geomHexa%T0, epsg%ECovaZETAZETA)+ &
                     matmul(geomHexa%TZETA, epsg%ECovaZETA))*zeta*zeta
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compGCovaMatrHexa
!
! Compute gradients for geometric matrix in covariant frame
!
! In  iNodeGeom        : first node of current geometric matrix
! In  jNodeGeom        : second node of current geometric matrix
! Out GCova0           : gradient for current geometric matrix (constant part)
! Out GCovaZETA        : gradient for current geometric matrix (ZETA part)
! Out GCovaZETAZETA    : gradient for current geometric matrix (ZETAZETA part)
!
! --------------------------------------------------------------------------------------------------
    subroutine compGCovaMatrHexa(iNodeGeom, jNodeGeom, GCova0, GCovaZETA, GCovaZETAZETA_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in)       :: iNodeGeom, jNodeGeom
        real(kind=8), intent(out) :: GCova0(SSH_SIZE_TENS)
        real(kind=8), intent(out) :: GCovaZETA(SSH_SIZE_TENS)
        real(kind=8), optional, intent(out) :: GCovaZETAZETA_(SSH_SIZE_TENS)
! - Local
        integer(kind=8) :: AD, EH, JM
        real(kind=8) :: aux13, aux23, aux33, auxz12, auxz13, auxz23
        real(kind=8) :: XI(3)
        real(kind=8) :: dN_dXsi(SSH_NBNODEG_HEXA), dN_dEta(SSH_NBNODEG_HEXA)
        real(kind=8) :: dN_dZeta(SSH_NBNODEG_HEXA)
        real(kind=8) :: GCovaZETAZETA(SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
        GCova0 = 0.d0
        GCovaZETA = 0.d0
        GCovaZETAZETA = 0.d0

! - Compute GCova0
        aux33 = 0.d0
        aux23 = 0.d0
        aux13 = 0.d0
        do AD = 1, 4
            XI = hexaQuadXIAD(AD, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux33 = aux33+dN_dZeta(iNodeGeom)*dN_dZeta(jNodeGeom)/4.d0
        end do
        do EH = 1, 4
            XI = hexaQuadXIEH(EH, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux23 = aux23+ &
                    (dN_dEta(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                     dN_dZeta(iNodeGeom)*dN_dEta(jNodeGeom))/4.d0
        end do
        do JM = 1, 4
            XI = hexaQuadXIJM(JM, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux13 = aux13+ &
                    (dN_dXsi(iNodeGeom)*dN_dZEta(jNodeGeom)+ &
                     dN_dZeta(iNodeGeom)*dN_dXsi(jNodeGeom))/4.d0
        end do
        GCova0(1) = hexaVectG1(iNodeGeom)*hexaVectG1(jNodeGeom)
        GCova0(2) = hexaVectG2(iNodeGeom)*hexaVectG2(jNodeGeom)
        GCova0(3) = aux33
        GCova0(4) = hexaVectG1(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                    hexaVectG2(iNodeGeom)*hexaVectG1(jNodeGeom)
        GCova0(5) = aux13
        GCova0(6) = aux23

! - Compute GCovaZETA
        auxz23 = 0.d0
        auxz13 = 0.d0
        do AD = 1, 4
            XI = hexaQuadXIEH(AD, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            auxz23 = auxz23+ &
                     XI(3)*(dN_dEta(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                            dN_dZeta(iNodeGeom)*dN_dEta(jNodeGeom))/4.d0
        end do
        do JM = 1, 4
            XI = hexaQuadXIJM(JM, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            auxz13 = auxz13+ &
                     XI(3)*(dN_dXsi(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                            dN_dZeta(iNodeGeom)*dN_dXsi(jNodeGeom))/4.d0
        end do
        auxz12 = hexaVectH2(iNodeGeom)*hexaVectG1(jNodeGeom)+ &
                 hexaVectH3(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                 hexaVectG2(iNodeGeom)*hexaVectH3(jNodeGeom)+ &
                 hexaVectG1(iNodeGeom)*hexaVectH2(jNodeGeom)
        GCovaZETA(1) = hexaVectH3(iNodeGeom)*hexaVectG1(jNodeGeom)+ &
                       hexaVectG1(iNodeGeom)*hexaVectH3(jNodeGeom)
        GCovaZETA(2) = hexaVectH2(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                       hexaVectG2(iNodeGeom)*hexaVectH2(jNodeGeom)
        GCovaZETA(3) = 0.d0
        GCovaZETA(4) = auxz12
        GCovaZETA(5) = auxz13
        GCovaZETA(6) = auxz23

! - Compute GCovaZETAZETA
        GCovaZETAZETA(1) = hexaVectH3(iNodeGeom)*hexaVectH3(jNodeGeom)
        GCovaZETAZETA(2) = hexaVectH2(iNodeGeom)*hexaVectH2(jNodeGeom)
        GCovaZETAZETA(3) = 0.d0
        GCovaZETAZETA(4) = hexaVectH2(iNodeGeom)*hexaVectH3(jNodeGeom)+ &
                           hexaVectH2(jNodeGeom)*hexaVectH3(iNodeGeom)
        GCovaZETAZETA(5) = 0.d0
        GCovaZETAZETA(6) = 0.d0

        if (present(GCovaZETAZETA_)) then
            GCovaZETAZETA_ = GCovaZETAZETA
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compGCovaSMatrHexa
!
! Compute gradients for geometric matrix of stabilization in covariant frame
!
! In  iNodeGeom        : first node of current geometric matrix
! In  jNodeGeom        : second node of current geometric matrix
! Out GCovaXI          : gradient for current geometric matrix (XI part)
! Out GCovaETA         : gradient for current geometric matrix (ETA part)
! Out GCovaETAZETA     : gradient for current geometric matrix (ETA/ZETA part)
! Out GCovaXIZETA      : gradient for current geometric matrix (XI/ZETA part)
!
! --------------------------------------------------------------------------------------------------
    subroutine compGCovaSMatrHexa(iNodeGeom, jNodeGeom, &
                                  GCovaXI, GCovaETA, &
                                  GCovaETAZETA, GCovaXIZETA)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in)       :: iNodeGeom, jNodeGeom
        real(kind=8), intent(out) :: GCovaXI(SSH_SIZE_TENS)
        real(kind=8), intent(out) :: GCovaETA(SSH_SIZE_TENS)
        real(kind=8), intent(out) :: GCovaETAZETA(SSH_SIZE_TENS)
        real(kind=8), intent(out) :: GCovaXIZETA(SSH_SIZE_TENS)
! - Local
        integer(kind=8) :: AD, EH, JM
        real(kind=8) :: aux13, aux23, aux33, aux11, aux12, aux22
        real(kind=8) :: XI(3)
        real(kind=8) :: dN_dXsi(SSH_NBNODEG_HEXA), dN_dEta(SSH_NBNODEG_HEXA)
        real(kind=8) :: dN_dZeta(SSH_NBNODEG_HEXA)
!   ------------------------------------------------------------------------------------------------
!
        GCovaXI = 0.d0
        GCovaETA = 0.d0
        GCovaETAZETA = 0.d0
        GCovaXIZETA = 0.d0

! - Compute GCovaXI
        aux33 = 0.d0
        aux23 = 0.d0
        do AD = 1, 4
            XI = hexaQuadXIAD(AD, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux33 = aux33+ &
                    XI(1)*(dN_dZeta(iNodeGeom)*dN_dZeta(jNodeGeom))/4.d0
        end do
        do EH = 1, 4
            XI = hexaQuadXIEH(EH, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux23 = aux23+ &
                    XI(1)*(dN_dEta(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                           dN_dZeta(iNodeGeom)*dN_dEta(jNodeGeom))/4.d0
        end do
        GCovaXI(1) = 0.d0
        GCovaXI(2) = hexaVectH1(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                     hexaVectG2(iNodeGeom)*hexaVectH1(jNodeGeom)
        GCovaXI(3) = aux33
        GCovaXI(4) = hexaVectH1(iNodeGeom)*hexaVectG1(jNodeGeom)+ &
                     hexaVectG1(iNodeGeom)*hexaVectH1(jNodeGeom)
        GCovaXI(5) = 0.d0
        GCovaXI(6) = aux23

! - Compute GCovaETA
        aux33 = 0.d0
        aux13 = 0.d0
        do AD = 1, 4
            XI = hexaQuadXIAD(AD, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux33 = aux33+ &
                    XI(2)*(dN_dZeta(iNodeGeom)*dN_dZeta(jNodeGeom))/4.d0
        end do
        do JM = 1, 4
            XI = hexaQuadXIJM(JM, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux13 = aux13+ &
                    XI(2)*(dN_dXsi(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                           dN_dZeta(iNodeGeom)*dN_dXsi(jNodeGeom))/4.d0
        end do
        GCovaETA(1) = hexaVectH1(iNodeGeom)*hexaVectG1(jNodeGeom)+ &
                      hexaVectG1(iNodeGeom)*hexaVectH1(jNodeGeom)
        GCovaETA(2) = 0.d0
        GCovaETA(3) = aux33
        GCovaETA(4) = hexaVectH1(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                      hexaVectG2(iNodeGeom)*hexaVectH1(jNodeGeom)
        GCovaETA(5) = aux13
        GCovaETA(6) = 0.d0

! - Compute GCovaETAZETA
        aux11 = hexaVectH4(iNodeGeom)*hexaVectG1(jNodeGeom)+ &
                hexaVectH3(iNodeGeom)*hexaVectH1(jNodeGeom)+ &
                hexaVectH1(iNodeGeom)*hexaVectH3(jNodeGeom)+ &
                hexaVectG1(iNodeGeom)*hexaVectH4(jNodeGeom)
        aux12 = hexaVectH4(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                hexaVectH2(iNodeGeom)*hexaVectH1(jNodeGeom)+ &
                hexaVectH1(iNodeGeom)*hexaVectH2(jNodeGeom)+ &
                hexaVectG2(iNodeGeom)*hexaVectH4(jNodeGeom)
        aux13 = 0.d0
        do JM = 1, 4
            XI = hexaQuadXIJM(JM, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux13 = aux13+ &
                    XI(2)*XI(3)*(dN_dXsi(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                                 dN_dXsi(jNodeGeom)*dN_dZeta(iNodeGeom))/4.d0
        end do
        GCovaETAZETA(1) = aux11
        GCovaETAZETA(2) = 0.d0
        GCovaETAZETA(3) = 0.d0
        GCovaETAZETA(4) = aux12
        GCovaETAZETA(5) = aux13
        GCovaETAZETA(6) = 0.d0

! - Compute GCovaXIZETA
        aux22 = hexaVectH4(iNodeGeom)*hexaVectG2(jNodeGeom)+ &
                hexaVectH2(iNodeGeom)*hexaVectH1(jNodeGeom)+ &
                hexaVectH1(iNodeGeom)*hexaVectH2(jNodeGeom)+ &
                hexaVectG2(iNodeGeom)*hexaVectH4(jNodeGeom)
        aux12 = hexaVectH4(iNodeGeom)*hexaVectG1(jNodeGeom)+ &
                hexaVectH3(iNodeGeom)*hexaVectH1(jNodeGeom)+ &
                hexaVectH1(iNodeGeom)*hexaVectH3(jNodeGeom)+ &
                hexaVectG1(iNodeGeom)*hexaVectH4(jNodeGeom)
        aux23 = 0.d0
        do EH = 1, 4
            XI = hexaQuadXIEH(EH, 1:3)
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
            aux23 = aux23+ &
                    XI(1)*XI(3)*(dN_dEta(iNodeGeom)*dN_dZeta(jNodeGeom)+ &
                                 dN_dEta(jNodeGeom)*dN_dZeta(iNodeGeom))/4.d0
        end do
        GCovaXIZETA(1) = 0.d0
        GCovaXIZETA(2) = aux22
        GCovaXIZETA(3) = 0.d0
        GCovaXIZETA(4) = aux12
        GCovaXIZETA(5) = 0.d0
        GCovaXIZETA(6) = aux23
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpslHexa
!
! Compute logarithmic strains in cartesian basis for HEXA cell
!
! In  epsg             : Green-Lagrange strains
! Out epsl             : logarithmic strains
! Out iretLog          : return code ==1 if logarithmic is not possible
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpslHexa(epsg, epsl, iretLog)
! - Parameters
        type(SSH_EPSG_HEXA), intent(in)  :: epsg
        type(SSH_EPSL_HEXA), intent(out) :: epsl
        integer(kind=8), intent(out)             :: iretLog
! - Local
        integer(kind=8), parameter :: nbvec = 3
        real(kind=8) :: epsc(6), epsl33(3, 3)
        integer(kind=8) :: i, j, k
!   ------------------------------------------------------------------------------------------------
!
        iretLog = 0
        epsl%vale = 0.d0
        epsl%eigenVale = 0.d0
        epsl%eigenVect = 0.d0
        epsl%logl = 0.d0

! - Vector form of strain tensor C=F'F=2E+I
        epsc(1) = 2.d0*epsg%vale(1)+1.d0
        epsc(2) = epsg%vale(4)
        epsc(3) = epsg%vale(5)
        epsc(4) = 2.d0*epsg%vale(2)+1.d0
        epsc(5) = epsg%vale(6)
        epsc(6) = 2.d0*epsg%vale(3)+1.d0

! - Eigen values/vectors of strain tensor
        call diagp3(epsc, epsl%eigenVect, epsl%eigenVale)
        do i = 1, nbvec
            if (epsl%eigenVale(i) .le. r8miem()) then
                iretLog = 1
                goto 999
            end if
            epsl%logl(i) = log(epsl%eigenVale(i))*0.5d0
        end do

! - Compute logarithmic strains
        epsl33 = 0.d0
        do i = 1, 3
            do j = 1, 3
                do k = 1, nbvec
                    epsl33(i, j) = epsl33(i, j)+ &
                                   epsl%logl(k)*epsl%eigenVect(i, k)*epsl%eigenVect(j, k)
                end do
            end do
        end do

! - Voigt form of strain tensor
        call tnsvec(3, 3, epsl33, epsl%vale, sqrt(2.d0))

999     continue
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Kinematic_Hexa_module
