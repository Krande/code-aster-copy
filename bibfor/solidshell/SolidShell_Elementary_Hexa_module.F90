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
public  :: compRigiMatrHexa, compSiefElgaHexa, compForcNodaHexa,&
           compRigiGeomHexaKpg,&
           compEpsgElgaHexa, compEpsiElgaHexa, compEpslElgaHexa,&
           compLoadHexa, compMassMatrHexa
! ==================================================================================================
private
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/SolidShell_type.h"
#include "asterfort/assert.h"
#include "asterfort/btsig.h"
#include "asterfort/jevecd.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
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
! --------------------------------------------------------------------------------------------------
!
! compForcNodaHexa
!
! Compute nodal forces for HEXA - FORC_NODA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  siefElga         : stresses at Gauss points
! Out forcNoda         : nodal forces
!
! --------------------------------------------------------------------------------------------------
subroutine compForcNodaHexa(elemProp, cellGeom,&
                            siefElga, forcNoda)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    real(kind=8), intent(in)        :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
    real(kind=8), intent(out)       :: forcNoda(SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
! - Local
    type(SSH_GEOM_HEXA) :: geomHexa
    type(SSH_KINE_HEXA) :: kineHexa
    real(kind=8) :: geomCurr(SSH_NBDOFG_HEXA), disp(SSH_NBDOF_HEXA)
    real(kind=8) :: zeta, poids, jacob
    integer :: nbIntePoint, kpg, jvCoor, jvWeight, jvCompor, jvDisp, iretc, iDof
    character(len=16) :: defoComp
!   ------------------------------------------------------------------------------------------------
!
    nbIntePoint = elemProp%elemInte%nbIntePoint
    jvCoor      = elemProp%elemInte%jvCoor
    jvWeight    = elemProp%elemInte%jvWeight

! - Prepare geometric quantities
    call initGeomCellHexa(cellGeom, geomHexa)
    if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! - Select configuration
    call tecach('ONO', 'PCOMPOR', 'L', iretc, iad = jvCompor)
    defoComp = 'PETIT'
    if (iretc .eq. 0) then
        defoComp = zk16(jvCompor-1+DEFO)
    endif

! - Update configuration
    if (defoComp .eq. 'PETIT') then
        geomCurr = cellGeom%geomInit
    else
        call tecach('ONO', 'PDEPLMR', 'L', iretc, iad = jvDisp)
        if (iretc .eq. 0) then
            do iDof = 1, SSH_NBDOF_HEXA
                disp(iDof) = zr(jvDisp-1+iDof)
            end do
            geomCurr = cellGeom%geomInit + disp(1:SSH_NBDOFG_HEXA)
        else
            call utmess('F', 'SOLIDSHELL1_4')
        endif
    endif

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

! ----- Product BT . sigma CUMULATED (see btsig subroutine)
        call btsig(elemProp%nbDof, SSH_SIZE_TENS, jacob,&
                   kineHexa%B, siefElga(1+SSH_SIZE_TENS*(kpg-1)), forcNoda)

    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRigiGeomHexaKpg
!
! Compute geometric matrix at current Gauss point for HEXA
!
! In  geomHexa         : geometric properties for HEXA cell
! In  cellGeom         : general geometric properties of cell
! In  matePara         : parameters of material
! Out matrRigi         : rigidity matrix
!
! --------------------------------------------------------------------------------------------------
subroutine compRigiGeomHexaKpg(geomHexa, zeta, sigm, matrGeom)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_GEOM_HEXA), intent(in) :: geomHexa
    real(kind=8), intent(in)        :: zeta, sigm(SSH_SIZE_TENS)
    real(kind=8), intent(out)       :: matrGeom(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
! - Local
    integer, parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
    integer :: iNodeGeom, jNodeGeom
    real(kind=8) :: const(SSH_SIZE_TENS)
    real(kind=8) :: GCart0(SSH_SIZE_TENS), GCartZETA(SSH_SIZE_TENS), GCartZETAZETA(SSH_SIZE_TENS)
    real(kind=8) :: GPinchZETA(SSH_SIZE_TENS), GPinchZZETA(SSH_SIZE_TENS), GPinchZZ(SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
    matrGeom = 0.d0

! - For "standard" nodes
    do iNodeGeom = 1, nbNodeGeom
        do jNodeGeom = 1, nbNodeGeom

! --------- Compute gradients for geometric matrix
            call compGCartMatrHexa(iNodeGeom, jNodeGeom, GCart0, GCartZETA)

! --------- Compute GCartZETAZETA
            GCartZETAZETA(1) = hexaVectH3(iNodeGeom)*hexaVectH3(jNodeGeom)
            GCartZETAZETA(2) = hexaVectH2(iNodeGeom)*hexaVectH2(jNodeGeom)
            GCartZETAZETA(3) = 0.d0
            GCartZETAZETA(4) = hexaVectH2(iNodeGeom)*hexaVectH3(jNodeGeom)+&
                               hexaVectH2(jNodeGeom)*hexaVectH3(iNodeGeom)
            GCartZETAZETA(5) = 0.d0
            GCartZETAZETA(6) = 0.d0

! --------- Compute matrix
            const = matmul(geomHexa%T, GCart0)+&
                    zeta*(matmul(geomHexa%T, GCartZETA) + &
                          matmul(geomHexa%TZETA, GCart0))+ &
                    zeta*zeta*(matmul(geomHexa%T, GCartZETAZETA) + &
                               matmul(geomHexa%TZETA, GCartZETA))
            matrGeom(3*(iNodeGeom-1)+1: 3*(iNodeGeom-1)+3, 3*(jNodeGeom-1)+1: 3*(jNodeGeom-1)+3) =&
                    sum(const*sigm)*matr3Iden
        enddo
    enddo

! - For "pinch" node
    GPinchZETA  = 0.d0
    GPinchZZETA = 0.d0
    GPinchZZ    = 0.d0
    do iNodeGeom = 1, nbNodeGeom
        GPinchZETA(3)  = -2.d0*hexaVectG3(iNodeGeom)
        GPinchZETA(5)  = -2.d0*hexaVectG1(iNodeGeom)
        GPinchZETA(6)  = -2.d0*hexaVectG2(iNodeGeom)
        GPinchZZETA(5) = -2.d0*hexaVectH3(iNodeGeom)
        GPinchZZETA(6) = -2.d0*hexaVectH2(iNodeGeom)
        const = zeta*matmul(geomHexa%T, GPinchZETA)+&
                zeta*zeta*(matmul(geomHexa%T, GPinchZZETA) + matmul(geomHexa%TZETA, GPinchZETA))
        matrGeom(25, 3*(iNodeGeom-1)+3) = sum(const*sigm)
        matrGeom(3*(iNodeGeom-1)+3, 25) = sum(const*sigm)
    enddo
    GPinchZZ(3) = 4.d0
    const = matmul(geomHexa%T, GPinchZZ)*zeta*zeta
    matrGeom(25, 25) = sum(const*sigm)
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsiElgaHexa
!
! Compute small strains for HEXA - EPSI_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  disp             : current displacements
! Out epsiElga         : small strains at Gauss points
!
! --------------------------------------------------------------------------------------------------
subroutine compEpsiElgaHexa(elemProp, cellGeom, disp, epsiElga)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
    real(kind=8), intent(out)       :: epsiElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
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
        epsiElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = epsi

    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsgElgaHexa
!
! Compute stresses for HEXA - EPSG_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  disp             : current displacements
! Out epsgElga         : Euler-Lagrange strains at Gauss points
!
! --------------------------------------------------------------------------------------------------
subroutine compEpsgElgaHexa(elemProp, cellGeom, disp, epsgElga)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
    real(kind=8), intent(out)       :: epsgElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! - Local
    type(SSH_GEOM_HEXA) :: geomHexa
    type(SSH_KINE_HEXA) :: kineHexa
    type(SSH_EPSG_HEXA) :: epsgHexa
    real(kind=8) :: zeta
    integer :: nbIntePoint, kpg, jvCoor
!   ------------------------------------------------------------------------------------------------
!
    nbIntePoint = elemProp%elemInte%nbIntePoint
    jvCoor      = elemProp%elemInte%jvCoor

! - Prepare geometric quantities
    call initGeomCellHexa(cellGeom, geomHexa)
    if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! - Flag for large strains
    kineHexa%lLarge = ASTER_TRUE

    do kpg = 1, nbIntePoint
        zeta  = zr(jvCoor-1+3*kpg)

! ----- Compute EAS B matrix in cartesian frame at current Gauss point
        call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! ----- Compute Green-Lagrange strains
        call compECovaMatrHexa(cellGeom, disp, epsgHexa)
        call compEpsgHexa(zeta, geomHexa, epsgHexa)
        epsgHexa%vale = epsgHexa%vale + kineHexa%BCartEAS * disp(25)
        if (SSH_DBG_KINE) call dbgObjEpsgHexa(epsgHexa)

        epsgElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = epsgHexa%vale(1:SSH_SIZE_TENS)

    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpslElgaHexa
!
! Compute stresses for HEXA - EPSL_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  disp             : current displacements
! Out epslElga         : logarithmic strains at Gauss points
!
! --------------------------------------------------------------------------------------------------
subroutine compEpslElgaHexa(elemProp, cellGeom,  disp, epslElga)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
    real(kind=8), intent(out)       :: epslElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! - Local
    type(SSH_GEOM_HEXA) :: geomHexa
    type(SSH_KINE_HEXA) :: kineHexa
    type(SSH_EPSG_HEXA) :: epsgHexa
    type(SSH_EPSL_HEXA) :: epslHexa
    real(kind=8) :: zeta
    integer :: cod, nbIntePoint, kpg, jvCoor
!   ------------------------------------------------------------------------------------------------
!
    nbIntePoint = elemProp%elemInte%nbIntePoint
    jvCoor      = elemProp%elemInte%jvCoor

! - Prepare geometric quantities
    call initGeomCellHexa(cellGeom, geomHexa)
    if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! - Flag for large strains
    kineHexa%lLarge = ASTER_TRUE

    do kpg = 1, nbIntePoint
        zeta  = zr(jvCoor-1+3*kpg)

! ----- Compute EAS B matrix in cartesian frame at current Gauss point
        call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! ----- Compute Green-Lagrange strains
        call compECovaMatrHexa(cellGeom, disp, epsgHexa)
        call compEpsgHexa(zeta, geomHexa, epsgHexa)
        epsgHexa%vale = epsgHexa%vale + kineHexa%BCartEAS * disp(25)
        if (SSH_DBG_KINE) call dbgObjEpsgHexa(epsgHexa)

! ----- Compute Logarithmic strains
        call compEpslHexa(epsgHexa, epslHexa, cod)
        if (SSH_DBG_KINE) call dbgObjEpslHexa(epslHexa)

        epslElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = epslHexa%vale

    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadHexa
!
! Compute loads for HEXA - CHAR_MECA_*
!
! In  cellGeom         : general geometric properties of cell
! In  option           : name of option to compute
! Out loadNoda         : nodal force from loads (Neumann)
!
! --------------------------------------------------------------------------------------------------
subroutine compLoadHexa(cellGeom, option, loadNoda)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    character(len=16), intent(in)   :: option
    real(kind=8), intent(out)       :: loadNoda(SSH_NBDOF_MAX)
! - Local
    integer :: jvPres
    real(kind=8) :: presSup, presInf, area
!   ------------------------------------------------------------------------------------------------
!
    loadNoda = 0.d0
    ASSERT(option .eq. 'CHAR_MECA_PRES_R')

! - Get input fields: for pressure, no node affected -> 0
    call jevecd('PPRESSR', jvPres, 0.d0)

! - Compute quantities in Ahmad frame for pinch quantities
    call compAhmadFrame(cellGeom, area)

! - Get pressure
    presSup = zr(jvPres-1+1)
    presInf = zr(jvPres-1+2)

! - Compute
    loadNoda(25) = loadNoda(25) +&
                   4.d0*(presInf-presSup)*area/3.d0
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMassMatrHexa
!
! Compute mass matrix for HEXA - MASS_MECA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  matePara         : parameters of material
! Out matrMass         : mass matrix
!
! --------------------------------------------------------------------------------------------------
subroutine compMassMatrHexa(elemProp, cellGeom, matePara, matrMass)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in) :: elemProp
    type(SSH_CELL_GEOM), intent(in) :: cellGeom
    type(SSH_MATE_PARA), intent(in) :: matePara
    real(kind=8), intent(out)       :: matrMass(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
! - Local
    integer, parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
    integer :: iNodeGeom, jNodeGeom
    real(kind=8) :: poids, jacob, XI(3)
    real(kind=8) :: rho(1), N(SSH_NBNODEG_HEXA), NPinch
    integer :: valeIret(1)
    character(len=4) :: inteFami
    integer :: nbIntePoint, kpg, jvCoor, jvWeight
    real(kind=8) :: matrMassPt(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
!
    matrMass    = 0.d0
    nbIntePoint = elemProp%elemInte%nbIntePoint
    inteFami    = elemProp%elemInte%inteFami
    jvCoor      = elemProp%elemInte%jvCoor
    jvWeight    = elemProp%elemInte%jvWeight

! - Loop on Gauss points
    do kpg = 1, nbIntePoint
        XI(1) = zr(jvCoor+3*(kpg-1)-1+1) 
        XI(2) = zr(jvCoor+3*(kpg-1)-1+2)
        XI(3) = zr(jvCoor+3*(kpg-1)-1+3)
        poids = zr(jvWeight-1+kpg)
        jacob = poids * cellGeom%detJac0

! ----- Get density
        call rcvalb(inteFami, kpg   , 1  , '+'        , matePara%jvMater,&
                    ' '     , 'ELAS', 0  , ' '        , [0.d0]          ,&
                    1       , 'RHO' , rho, valeIret(1), 1)

! ----- Construction of shape functions
        N  = hexaVectS1 + &
             XI(1)*hexaVectG1+XI(2)*hexaVectG2+XI(3)*hexaVectG3+&
             XI(1)*XI(2)*hexaVectH1+&
             XI(3)*XI(2)*hexaVectH2+&
             XI(1)*XI(3)*hexaVectH3+&
             XI(1)*XI(2)*XI(3)*hexaVectH4
        NPinch = 1.d0 - XI(3)*XI(3)

! ----- Compute matrix on volumic nodes
        matrMassPt = 0.d0
        do iNodeGeom = 1, nbNodeGeom
            do jNodeGeom = 1, nbNodeGeom
                matrMassPt(3*(iNodeGeom-1)+1:3*(iNodeGeom-1)+3,&
                           3*(jNodeGeom-1)+1:3*(jNodeGeom-1)+3) = &
                    rho(1) * jacob * N(iNodeGeom) * N(jNodeGeom) * matr3Iden
            end do
        end do

! ----- Compute matrix on pinch node
        do iNodeGeom = 1, nbNodeGeom
            matrMassPt(SSH_NBDOF_HEXA, 3*(iNodeGeom-1)+3) = rho(1) * jacob*N(iNodeGeom) * NPinch
            matrMassPt(3*(iNodeGeom-1)+3, SSH_NBDOF_HEXA) = rho(1) * jacob*N(iNodeGeom) * NPinch
        enddo
        matrMassPt(SSH_NBDOF_HEXA, SSH_NBDOF_HEXA) = rho(1) * poids * NPinch * NPinch

! ----- Update matrix
        matrMass = matrMass + matrMassPt

    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
!
end module SolidShell_Elementary_Hexa_module
