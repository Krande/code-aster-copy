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
subroutine irmano(meshNameZ, &
                  nbCellSelect, cellSelect, &
                  nodeFlag)
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
!
    character(len=*), intent(in) :: meshNameZ
    integer(kind=8), intent(in) :: nbCellSelect
    integer(kind=8), pointer :: cellSelect(:)
    aster_logical, pointer :: nodeFlag(:)
!
! --------------------------------------------------------------------------------------------------
!
! Print results
!
! Select nodes from cells
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: cellNbNode
    character(len=8) :: meshName
    integer(kind=8) :: iCell, cellNume, iNode, nodeFirst, nodeNume
    integer(kind=8), pointer :: conxLong(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    meshName = meshNameZ
!
! - Access to connectivity
!
    call jeveuo(meshName//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(meshName//'.CONNEX', 'LONCUM'), 'L', vi=conxLong)
!
! - Loop on cells
!
    do iCell = 1, nbCellSelect
        cellNume = cellSelect(iCell)
        nodeFirst = conxLong(cellNume)
        cellNbNode = conxLong(cellNume+1)-nodeFirst
        do iNode = 1, cellNbNode
            nodeNume = connex(nodeFirst-1+iNode)
            if (.not. nodeFlag(nodeNume)) then
                nodeFlag(nodeNume) = ASTER_TRUE
            end if
        end do
    end do
!
    call jedema()
end subroutine
