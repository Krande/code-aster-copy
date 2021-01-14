! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! Module for elementary computations of HEXA cells for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Elementary_Hexa_module
! ==================================================================================================
use SolidShell_type
use SolidShell_Utilities_module
use SolidShell_Debug_module
use SolidShell_Geometry_Hexa_module
use SolidShell_Kinematic_Hexa_module
use SolidShell_Stabilization_Hexa_module
! ==================================================================================================
implicit none
! ==================================================================================================
public  :: compRigiMatrHexa, compSiefElgaHexa
! ==================================================================================================
private
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/SolidShell_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! compRigiMatrHexa
!
! Compute rigidity matrix for HEXA - RIGI_MECA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  matePara         : parameters of material
! Out matrRigi         : rigidity matrix
!
! --------------------------------------------------------------------------------------------------
subroutine compRigiMatrHexa(elemProp, cellGeom, matePara, matrRigi)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    type(SSH_MATE_PARA), intent(in) :: matePara
    real(kind=8), intent(out)       :: matrRigi(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
! - Local
    type(SSH_GEOM_HEXA) :: geomHexa
    type(SSH_KINE_HEXA) :: kineHexa
    type(SSH_STAB_HEXA) :: stabHexa
    real(kind=8) :: geomCurr(SSH_NBDOFG_HEXA)
    real(kind=8) :: zeta, poids, jacob, Ueff
    integer :: nbIntePoint, kpg, jvCoor, jvWeight
    real(kind=8) :: tBDB(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
!
    nbIntePoint = elemProp%elemInte%nbIntePoint
    jvCoor      = elemProp%elemInte%jvCoor
    jvWeight    = elemProp%elemInte%jvWeight

! - Prepare geometric quantities
    call initGeomCellHexa(cellGeom, geomHexa)
    if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! - Update configuration
    geomCurr = cellGeom%geomInit

! - Compute gradient matrix in covariant basis
    call compBCovaMatrHexa(geomCurr, kineHexa)

! - Compute gradient matrix in cartesian frame
    call compBCartMatrHexa(geomHexa, kineHexa)

    if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_ = ASTER_TRUE)

! - Loop on Gauss points
    do kpg = 1, nbIntePoint
        zeta  = zr(jvCoor-1+3*kpg)
        poids = zr(jvWeight-1+kpg)
        jacob = poids * cellGeom%detJac0

! ----- Compute EAS B matrix in cartesian frame at current Gauss point
        call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! ----- Compute B matrix
        call compBMatrHexa(zeta, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_ = ASTER_TRUE)

! ----- Compute product tBSB
        call prodBTDB(matePara%elemHookeMatrix, SSH_SIZE_TENS, elemProp%nbDof, kineHexa%B, tBDB)

! ----- Update matrix
        matrRigi = matrRigi +  jacob * tBDB
        if (SSH_DBG_ELEM) call dbgMatrRigiKpg(kpg, zeta, poids, jacob, matrRigi)

    end do

! - Effective shear modulus for stabilization (elasticity)
    Ueff = matePara%elemHookeMatrix(5, 5)

! - Compute stabilization matrix (material part, elasticity)
    call compStabMatrMateHexa(geomHexa, kineHexa, Ueff, stabHexa)
    if (SSH_DBG_STAB) call dbgObjStabHexa(Ueff, stabHexa)

! - Compute matrix
    matrRigi(1:SSH_NBDOFG_HEXA, 1:SSH_NBDOFG_HEXA) = &
        matrRigi(1:SSH_NBDOFG_HEXA, 1:SSH_NBDOFG_HEXA) + stabHexa%matrStabMate
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSiefElgaHexa
!
! Compute stresses for HEXA - SIEF_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  matePara         : parameters of material
! In  disp             : current displacements
! Out siefElga         : stresses at Gauss points
!
! --------------------------------------------------------------------------------------------------
subroutine compSiefElgaHexa(elemProp, cellGeom, matePara, disp,&
                            siefElga)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    type(SSH_MATE_PARA), intent(in) :: matePara
    real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
    real(kind=8), intent(out)       :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! - Local
    type(SSH_GEOM_HEXA) :: geomHexa
    type(SSH_KINE_HEXA) :: kineHexa
    real(kind=8) :: geomCurr(SSH_NBDOFG_HEXA)
    real(kind=8) :: zeta, epsi(SSH_SIZE_TENS)
    integer :: nbIntePoint, kpg, jvCoor
!   ------------------------------------------------------------------------------------------------
!
    nbIntePoint = elemProp%elemInte%nbIntePoint
    jvCoor      = elemProp%elemInte%jvCoor

! - Prepare geometric quantities
    call initGeomCellHexa(cellGeom, geomHexa)
    if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! - Update configuration
    geomCurr = cellGeom%geomInit

! - Compute gradient matrix in covariant basis
    call compBCovaMatrHexa(geomCurr, kineHexa)

! - Compute gradient matrix in cartesian frame
    call compBCartMatrHexa(geomHexa, kineHexa)
    if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_ = ASTER_TRUE)

! - Loop on Gauss points
    do kpg = 1, nbIntePoint
        zeta  = zr(jvCoor-1+3*kpg)

! ----- Compute EAS B matrix in cartesian frame at current Gauss point
        call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! ----- Compute B matrix
        call compBMatrHexa(zeta, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_ = ASTER_TRUE)

! ----- Compute small strains
        call compEpsiHexa(kineHexa, disp, epsi)

! ----- Compute stresses
        siefElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = matmul(matePara%elemHookeMatrix, epsi)

    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
!
end module SolidShell_Elementary_Hexa_module
