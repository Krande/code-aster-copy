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
subroutine irmama(meshNameZ, &
                  nbCell, cellName, &
                  nbGrCell, grCellName, &
                  nbCellSelect, cellFlag, lfichUniq)
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*), intent(in) :: meshNameZ
    integer(kind=8), intent(in) :: nbCell
    character(len=8), pointer :: cellName(:)
    integer(kind=8), intent(in) :: nbGrCell
    character(len=24), pointer :: grCellName(:)
    integer(kind=8), intent(out) :: nbCellSelect
    aster_logical, pointer :: cellFlag(:)
    aster_logical, intent(in) :: lfichUniq
!
! --------------------------------------------------------------------------------------------------
!
! Print results
!
! Select cells from user
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: meshName
    character(len=11) :: vecGrpName
    integer(kind=8) :: iCell, cellNume, iGrCell, iret, grCellNbCell
    integer(kind=8), pointer :: listCell(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    meshName = meshNameZ
    nbCellSelect = 0
!
! - Select cells by name
!
    if (nbCell .ne. 0) then
        if (lfichUniq) then
            call utmess('F', 'MED3_4')
        end if
        do iCell = 1, nbCell
            cellNume = char8_to_int(cellName(iCell))
            if (cellNume .eq. 0) then
                call utmess('A', 'RESULT3_9', sk=cellName(iCell))
                cellName(iCell) = ' '
            else
                if (.not. cellFlag(cellNume)) then
                    cellFlag(cellNume) = ASTER_TRUE
                    nbCellSelect = nbCellSelect+1
                end if
            end if
        end do
    end if
!
! - Select cells in groups of cells
!
    if (nbGrCell .ne. 0) then
        vecGrpName = '.GROUPEMA'
        if (lfichUniq) vecGrpName = '.PAR_GRPMAI'
        do iGrCell = 1, nbGrCell
            call jeexin(jexnom(meshName//vecGrpName, grCellName(iGrCell)), iret)
            if (iret .eq. 0) then
                call utmess('A', 'RESULT3_10', sk=grCellName(iGrCell))
                grCellName(iGrCell) = ' '
            else
                call jeexin(jexnom(meshName//'.GROUPEMA', grCellName(iGrCell)), iret)
                if (iret .ne. 0) then
                    call jelira(jexnom(meshName//'.GROUPEMA', grCellName(iGrCell)), &
                                'LONMAX', grCellNbCell)
                    if (grCellNbCell .eq. 0) then
                        call utmess('A', 'RESULT3_11', sk=grCellName(iGrCell))
                        grCellName(iGrCell) = ' '
                    else
                        call jeveuo(jexnom(meshName//'.GROUPEMA', grCellName(iGrCell)), &
                                    'L', vi=listCell)
                        do iCell = 1, grCellNbCell
                            cellNume = listCell(iCell)
                            if (cellNume .ne. 0) then
                                if (.not. cellFlag(cellNume)) then
                                    cellFlag(cellNume) = ASTER_TRUE
                                    nbCellSelect = nbCellSelect+1
                                end if
                            end if
                        end do
                    end if
                end if
            end if
        end do
    end if
!
    call jedema()
end subroutine
