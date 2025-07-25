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
subroutine resuSelectNode(meshName, meshNodeNb, &
                          nodeUserNb, nodeUserNume, &
                          nodeName, nodeNume, &
                          nodeNb)
!
    implicit none
!
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: meshName
    integer(kind=8), intent(in) :: nodeUserNb, meshNodeNb, nodeUserNume(*)
    character(len=8), pointer :: nodeName(:)
    integer(kind=8), pointer :: nodeNume(:)
    integer(kind=8), intent(out) :: nodeNb
!
! --------------------------------------------------------------------------------------------------
!
! Result management
!
! Select list of nodes
!
! --------------------------------------------------------------------------------------------------
!
! In  meshName         : name of mesh
! In  meshNodeNb       : number of nodes require in mesh
! In  nodeUserNb       : number of nodes require by user
! In  nodeUserNume     : list of index of nodes require by user
! Ptr nodeName         : pointer to the name of nodes selected
! Ptr nodeNume         : pointer to the index of nodes selected
! Out nodeNb           : number of nodes selected
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iNode
!
! --------------------------------------------------------------------------------------------------
!
    nodeNb = 0
    if (nodeUserNb .eq. 0) then
        do iNode = 1, meshNodeNb
            nodeName(iNode) = int_to_char8(iNode)
            nodeNume(iNode) = iNode
            nodeNb = meshNodeNb
        end do
    else
        do iNode = 1, nodeUserNb
            nodeNume(iNode) = nodeUserNume(iNode)
            nodeName(nodeNume(iNode)) = int_to_char8(iNode)
        end do
        nodeNb = nodeUserNb
    end if
!
end subroutine
