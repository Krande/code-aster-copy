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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine vetyma(mesh, ndim, loadType, listCell, nbCell)
!
    use mesh_module, only: getPropertiesOfCell
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbCell
    character(len=24), intent(in) :: listCell
    character(len=16), intent(in) :: loadType
    integer(kind=8), intent(in) :: ndim
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Check element type
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh      : name of mesh
! In  ndim      : space dimension
! In  loadType  : type of load
! In  listCell  : list of elements read
! In  nbCell    : number of elements read
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: valk(2)
    character(len=8) :: cellName, cellTypeName
    character(len=4) :: cellTopo, topo_2d, topo_3d, topoRequired
    integer(kind=8) :: nerr
    integer(kind=8) :: cellNume, cellOrder, cellTypeNume
    integer(kind=8) :: iCell, orderMini
    integer(kind=8), pointer :: listCellNume(:) => null()
    integer(kind=8), pointer :: meshTypmail(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3 .or. ndim .eq. 23)
    if (nbCell .eq. 0) goto 99
!
! - Access to mesh
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypmail)
    call jeveuo(listCell, 'L', vi=listCellNume)
!
! - Type of elements
!
    orderMini = -1
    if (loadType .eq. 'FLUX_REP' .or. loadType .eq. 'PRES_REP' .or. &
        loadType .eq. 'ECHANGE' .or. loadType .eq. 'FORCE_FACE' .or. &
        loadType .eq. 'VITE_FACE' .or. &
        loadType .eq. 'FORCE_CONTOUR' .or. loadType .eq. 'EFFE_FOND' .or. &
        loadType .eq. 'ONDE_PLAN') then
        topo_2d = 'LINE'
        topo_3d = 'SURF'

    else if (loadType .eq. 'SOURCE' .or. loadType .eq. 'FORCE_INTERNE') then
        topo_2d = 'SURF'
        topo_3d = 'VOLU'

    else if (loadType .eq. 'FORCE_TUYAU') then
        topo_2d = 'None'
        topo_3d = 'LINE'
        orderMini = 2

    else
        goto 99
    end if
!
! - Select the right set of topological dimension of cells
!
    if (ndim .eq. 2) then
        topoRequired = topo_2d
    else if (ndim .eq. 3) then
        topoRequired = topo_3d
    else if (ndim .eq. 23) then
! ----- We determine from the first element
        cellNume = listCellNume(1)
        cellName = int_to_char8(cellNume)
        cellTypeNume = meshTypmail(cellNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
        call getPropertiesOfCell(cellTypeName, topoRequired)
    end if
!
    nerr = 0
    do iCell = 1, nbCell

! ----- Current cell
        cellNume = listCellNume(iCell)
        cellName = int_to_char8(cellNume)
        cellTypeNume = meshTypmail(cellNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)

! ----- Determine properties of cell
        call getPropertiesOfCell(cellTypeName, cellTopo, cellOrder)

! ----- Check consistency of topological shape
        if (topoRequired .eq. 'LINE') then
            if (cellTopo .ne. 'LINE') then
                nerr = nerr+1
                valk(1) = cellName
                valk(2) = loadType
                call utmess('A', 'CHARGES2_86', nk=2, valk=valk)
            end if
        else if (topoRequired .eq. 'SURF') then
            if (cellTopo .ne. 'SURF') then
                nerr = nerr+1
                valk(1) = cellName
                valk(2) = loadType
                call utmess('A', 'CHARGES2_87', nk=2, valk=valk)
            end if
        else if (topoRequired .eq. 'VOLU') then
            if (cellTopo .ne. 'VOLU') then
                nerr = nerr+1
                valk(1) = cellName
                valk(2) = loadType
                call utmess('A', 'CHARGES2_88', nk=2, valk=valk)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Check consistency of interpolation order
        if (orderMini .ge. 0) then
            if (cellOrder .lt. orderMini) then
                valk(1) = cellName
                valk(2) = loadType
                call utmess('A', 'CHARGES2_85', nk=2, valk=valk)
            end if
        end if

    end do
!
    if (nbCell .eq. nerr) then
        call utmess('A', 'CHARGES2_89', sk=loadType)
    end if
!
99  continue
    call jedema()
end subroutine
