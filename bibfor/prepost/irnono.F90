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
subroutine irnono(meshNameZ, &
                  nbNode, nodeName, &
                  nbGrNode, grNodeName, &
                  nbNodeSelect, nodeFlag, lfichUniq)
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
    integer(kind=8), intent(in) :: nbNode
    character(len=8), pointer :: nodeName(:)
    integer(kind=8), intent(in) :: nbGrNode
    character(len=24), pointer :: grNodeName(:)
    integer(kind=8), intent(out) :: nbNodeSelect
    aster_logical, pointer :: nodeFlag(:)
    aster_logical, intent(in) :: lfichUniq
!
! --------------------------------------------------------------------------------------------------
!
! Print results
!
! Select nodes from user
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: meshName
    character(len=11) :: vecGrpName
    integer(kind=8) :: iNode, nodeNume, iGrNode, iret, grNodeNbNode
    integer(kind=8), pointer :: listNode(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    meshName = meshNameZ
    nbNodeSelect = 0
!
! - Select nodes by name
!
    if (nbNode .ne. 0) then
        if (lfichUniq) then
            call utmess('F', 'MED3_4')
        end if
        do iNode = 1, nbNode
            nodeNume = char8_to_int(nodeName(iNode))
            if (nodeNume .eq. 0) then
                call utmess('A', 'RESULT3_6', sk=nodeName(iNode))
                nodeName(iNode) = ' '
            else
                if (.not. nodeFlag(nodeNume)) then
                    nodeFlag(nodeNume) = ASTER_TRUE
                    nbNodeSelect = nbNodeSelect+1
                end if
            end if
        end do
    end if
!
! - Select nodes in groups of nodes
!
    if (nbGrNode .ne. 0) then
        vecGrpName = '.GROUPENO'
        if (lfichUniq) vecGrpName = '.PAR_GRPNOE'
        do iGrNode = 1, nbGrNode
            call jeexin(jexnom(meshName//vecGrpName, grNodeName(iGrNode)), iret)
            if (iret .eq. 0) then
                call utmess('A', 'RESULT3_7', sk=grNodeName(iGrNode))
                grNodeName(iGrNode) = ' '
            else
                call jeexin(jexnom(meshName//'.GROUPENO', grNodeName(iGrNode)), iret)
                if (iret .ne. 0) then
                    call jelira(jexnom(meshName//'.GROUPENO', grNodeName(iGrNode)), &
                                'LONMAX', grNodeNbNode)
                    if (grNodeNbNode .eq. 0) then
                        call utmess('A', 'RESULT3_8', sk=grNodeName(iGrNode))
                        grNodeName(iGrNode) = ' '
                    else
                        call jeveuo(jexnom(meshName//'.GROUPENO', grNodeName(iGrNode)), &
                                    'L', vi=listNode)
                        do iNode = 1, grNodeNbNode
                            nodeNume = listNode(iNode)
                            if (.not. nodeFlag(nodeNume)) then
                                nodeFlag(nodeNume) = ASTER_TRUE
                                nbNodeSelect = nbNodeSelect+1
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
