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
module model_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: getFSICell, getFluidCell, getPlateCell
    private :: getAccess
! ==================================================================================================
    private
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexatr.h"
#include "asterfort/lteatt.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! getAccess
!
! Get access to mesh and model
!
! In  model            : name of model
! Out nbCell           : number of cells in mesh
! Ptr modelCells       : link between cell and finite element
! Out hasFE            : flag for finite element
!
! --------------------------------------------------------------------------------------------------
    subroutine getAccess(modelz, nbCell, modelCells, hasFE)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelz
        integer, intent(out) :: nbCell
        integer, pointer :: modelCells(:)
        aster_logical, intent(out) :: hasFE
! ----- Local
        integer :: iret
        character(len=8) :: model, mesh
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        nbCell = 0
        modelCells => null()
        hasFE = ASTER_FALSE
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)
        call jeexin(model//'.MAILLE', iret)
        if (iret .ne. 0) then
            hasFE = ASTER_TRUE
            call jeveuo(model//'.MAILLE', 'L', vi=modelCells)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getFSICell
!
! Get list of cells with FSI model
!
! In  model            : name of model
! Out nbCellFSI        : number of FSI cells
! OUT cellFSI          : for each cell in mesh a flag for FSI cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getFSICell(modelz, nbCellFSI, cellFSI)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelz
        integer, intent(out) :: nbCellFSI
        integer, pointer :: cellFSI(:)
! ----- Local
        integer :: iCell, nbCell, iCellFSI
        integer :: cellTypeNume
        integer, pointer :: modelCells(:) => null()
        character(len=16) :: cellTypeName
        aster_logical, pointer :: isCellFSI(:) => null()
        aster_logical :: hasFE
!   ------------------------------------------------------------------------------------------------
!
        nbCellFSI = 0
        call getAccess(modelz, nbCell, modelCells, hasFE)
        if (hasFE) then
            AS_ALLOCATE(vl=isCellFSI, size=nbCell)
            do iCell = 1, nbCell
                cellTypeNume = modelCells(iCell)
                if (cellTypeNume .ne. 0) then
                    call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
                    if (lteatt('FSI', 'OUI', typel=cellTypeName)) then
                        nbCellFSI = nbCellFSI+1
                        isCellFSI(iCell) = ASTER_TRUE
                    end if
                end if
            end do
            if (nbCellFSI .ne. 0) then
                iCellFSI = 0
                AS_ALLOCATE(vi=cellFSI, size=nbCellFSI)
                do iCell = 1, nbCell
                    if (isCellFSI(iCell)) then
                        iCellFSI = iCellFSI+1
                        cellFSI(iCellFSI) = iCell
                    end if
                end do
                ASSERT(iCellFSI .eq. nbCellFSI)
            end if
            AS_DEALLOCATE(vl=isCellFSI)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getFluidCell
!
! Get list of cells with fluid model
!
! In  model            : name of model
! Out nbCellFluid      : number of fluid cells
! OUT cellFluid        : for each cell in mesh a flag for fluid cells
!
! --------------------------------------------------------------------------------------------------
    subroutine getFluidCell(modelz, nbCellFluid, cellFluid)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelz
        integer, intent(out) :: nbCellFluid
        integer, pointer :: cellFluid(:)
! ----- Local
        integer :: iCell, nbCell, iCellFluid
        integer :: cellTypeNume
        integer, pointer :: modelCells(:) => null()
        character(len=16) :: cellTypeName
        aster_logical, pointer :: isCellFluid(:) => null()
        aster_logical :: hasFE
!   ------------------------------------------------------------------------------------------------
!
        nbCellFluid = 0
        call getAccess(modelz, nbCell, modelCells, hasFE)
        if (hasFE) then
            AS_ALLOCATE(vl=isCellFluid, size=nbCell)
            do iCell = 1, nbCell
                cellTypeNume = modelCells(iCell)
                if (cellTypeNume .ne. 0) then
                    call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
                    if (lteatt('FLUIDE', 'OUI', typel=cellTypeName)) then
                        nbCellFluid = nbCellFluid+1
                        isCellFluid(iCell) = ASTER_TRUE
                    end if
                end if
            end do
            iCellFluid = 0
            AS_ALLOCATE(vi=cellFluid, size=nbCellFluid)
            do iCell = 1, nbCell
                if (isCellFluid(iCell)) then
                    iCellFluid = iCellFluid+1
                    cellFluid(iCellFluid) = iCell
                end if
            end do
            ASSERT(iCellFluid .eq. nbCellFluid)
            AS_DEALLOCATE(vl=isCellFluid)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getPlateCell
!
! Get list of cells with plate model
!
! In  model            : name of model
! Out nbCellPlate      : number of plate cells
! Ptr cellPlate        : pointer to cells with plate model
!
! --------------------------------------------------------------------------------------------------
    subroutine getPlateCell(modelz, nbCellPlate, cellPlate)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelz
        integer, intent(out) :: nbCellPlate
        integer, pointer :: cellPlate(:)
! ----- Local
        integer :: iCell, nbCell, iCellPlate
        integer :: cellTypeNume
        integer, pointer :: modelCells(:) => null()
        character(len=16) :: cellTypeName
        aster_logical, pointer :: isCellPlate(:) => null()
        aster_logical :: hasFE
!   ------------------------------------------------------------------------------------------------
!
        nbCellPlate = 0
        call getAccess(modelz, nbCell, modelCells, hasFE)
        if (hasFE) then
            AS_ALLOCATE(vl=isCellPlate, size=nbCell)
            do iCell = 1, nbCell
                cellTypeNume = modelCells(iCell)
                if (cellTypeNume .ne. 0) then
                    call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
                    if (lteatt('PLAQUE', 'OUI', typel=cellTypeName)) then
                        nbCellPlate = nbCellPlate+1
                        isCellPlate(iCell) = ASTER_TRUE
                    end if
                end if
            end do
            if (nbCellPlate .ne. 0) then
                iCellPlate = 0
                AS_ALLOCATE(vi=cellPlate, size=nbCellPlate)
                do iCell = 1, nbCell
                    if (isCellPlate(iCell)) then
                        iCellPlate = iCellPlate+1
                        cellPlate(iCellPlate) = iCell
                    end if
                end do
                ASSERT(iCellPlate .eq. nbCellPlate)
            end if
            AS_DEALLOCATE(vl=isCellPlate)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module model_module
