! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine te0425(option, nomte)
!
    use contact_module
    use contact_type
    use ContactPairing_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "contact_module.h"
#include "Contact_type.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!  Compute Geometric Gap for Mortar methods
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    integer :: jvGeom, jvGapVale, jvGapStatus
    integer :: iNodeSlav
    real(kind=8) :: normSlav(3), coorProjGlob(3), gap
    type(Contact_CellGeom) :: contCellGeom
    type(Cell_Geom) :: cellMastLine
    type(Contact_ProjAlgoPara) :: projAlgoPara
    type(Contact_ProjPara) :: projOnMast
    aster_logical :: pointIsInside
!
! --------------------------------------------------------------------------------------------------
!

! - Get input and output fields
    call jevech('PGEOMCR', 'L', jvGeom)
    call jevech('PVECGAP', 'E', jvGapVale)
    call jevech('PVEIGAP', 'E', jvGapStatus)

! - Set parameters of Newton algorithm for projection
    call setProjAlgoPara(PROJ_ALGO_GAPI, debug, projAlgoPara)

! - Get geometry of contact cell
    call getContactCellGeom(nomte, jvGeom, contCellGeom)

! - No values on master side
    zr(jvGapVale-1+1:jvGapVale-1+contCellGeom%slav%nbNode) = 0.d0
    zr(jvGapStatus-1+1:jvGapStatus-1+contCellGeom%slav%nbNode) = 0.d0
!
    if (contCellGeom%slav%cellCode .ne. "POI1") then

! ----- Linearize master cell
        call linearizeCell(contCellGeom%mast, cellMastLine)

! ----- Prepare object for projection on master cell
        projOnMast%modelDime = contCellGeom%cellDime
        projOnMast%geomTarget = contCellGeom%mast
        projOnMast%geomTargetLine = cellMastLine
        ASSERT(projAlgoPara%withPrepLine)

! ----- Projection of all slave nodes
        do iNodeSlav = 1, contCellGeom%slav%nbNode

! --------- Compute external norm to slave cell at current slave node
            call calcNormExte(contCellGeom%slav, contCellGeom%cellDime, iNodeSlav, normSlav)

! --------- Projection of slave node on master cell
            projOnMast%pointCoor = contCellGeom%slav%coorNodeGlob(:, iNodeSlav)
            projOnMast%projVect = normSlav
            call projPointVector(projAlgoPara, projOnMast, pointIsInside)
            if (.not. pointIsInside) then
                cycle
            end if
            ASSERT(projOnMast%errorCode == 0)

! --------- Return in real master space
            call fromParaToGLob(contCellGeom%mast, contCellGeom%cellDime, &
                                projOnMast%ksi, coorProjGlob)

! --------- Compute raytracing gap
            gap = gapEval(contCellGeom%slav%coorNodeGlob(:, iNodeSlav), coorProjGlob, normSlav)

! --------- Save distance and status
            zr(jvGapVale-1+iNodeSlav) = gap
            zr(jvGapStatus-1+iNodeSlav) = 1.d0

        end do
    end if
!
end subroutine
