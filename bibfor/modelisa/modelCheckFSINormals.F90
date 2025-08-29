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
subroutine modelCheckFSINormals(model)
!
    use mesh_module, only: getPropertiesOfListOfCells, getSkinCellSupport, checkNormalOnSkinCell
    use model_module, only: getFSICell, getFluidCell
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_MODELE
!
! Check normals for FSI elements
!
! --------------------------------------------------------------------------------------------------
!
! In  model           : name of the model
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: mesh
    integer(kind=8) :: nbCell, modelDime
    integer(kind=8) :: nbCellFSI, nbCellFluid
    integer(kind=8), pointer :: cellFSI(:) => null()
    integer(kind=8), pointer :: cellFluid(:) => null()
    aster_logical :: lCellSurf, lCellLine
    integer(kind=8), pointer :: cellFSINbNode(:) => null()
    integer(kind=8), pointer :: cellFSINodeIndx(:) => null()
    integer(kind=8), pointer :: cellFSISupp(:) => null()
    aster_logical :: lMisoriented
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)

    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)
    call dismoi('DIM_GEOM', model, 'MODELE', repi=modelDime)
!
! - Get list of cells with FSI model
!
    call getFSICell(model, nbCellFSI, cellFSI)

    if (nbCellFSI .gt. 0) then
        if (niv .ge. 2) then
            call utmess('I', 'MODELE1_80', si=nbCellFSI)
        end if
! ----- Get list of cells with fluid model
        call getFluidCell(model, nbCellFluid, cellFluid)
! ----- Get properties of FSI cells
        AS_ALLOCATE(vi=cellFSINbNode, size=nbCellFSI)
        AS_ALLOCATE(vi=cellFSINodeIndx, size=nbCellFSI)
        call getPropertiesOfListOfCells(mesh, &
                                        nbCellFSI, cellFSI, &
                                        cellFSINbNode, cellFSINodeIndx, &
                                        lCellSurf, lCellLine)
        ASSERT(lCellLine .or. lCellSurf)
        if (lCellLine) then
            ASSERT(.not. lCellSurf)
            if (niv .ge. 2) then
                call utmess('I', 'MODELE1_81')
            end if
        end if
        if (lCellSurf) then
            ASSERT(.not. lCellLine)
            if (niv .ge. 2) then
                call utmess('I', 'MODELE1_82')
            end if
        end if

! ----- Get "volumic" cells support of skin cells
        AS_ALLOCATE(vi=cellFSISupp, size=nbCellFSI)
        call getSkinCellSupport(mesh, &
                                nbCellFSI, cellFSI, &
                                lCellSurf, lCellLine, &
                                cellFSISupp, &
                                nbCellFluid, cellFluid)

! ----- Check normals
        call checkNormalOnSkinCell(mesh, modelDime, &
                                   nbCellFSI, cellFSI, &
                                   cellFSINbNode, cellFSINodeIndx, &
                                   cellFSISupp, lMisoriented)
        if (lMisoriented) then
            call utmess('F', 'FLUID1_4')
        end if
!
        AS_DEALLOCATE(vi=cellFSINbNode)
        AS_DEALLOCATE(vi=cellFSINodeIndx)
        AS_DEALLOCATE(vi=cellFSISupp)
    end if

!
    AS_DEALLOCATE(vi=cellFSI)
    AS_DEALLOCATE(vi=cellFluid)
!
end subroutine
