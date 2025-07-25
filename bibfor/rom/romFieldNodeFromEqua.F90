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
subroutine romFieldNodeFromEqua(field, nbNodeMesh, listNode)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_Field), intent(in) :: field
    integer(kind=8), intent(in) :: nbNodeMesh
    integer(kind=8), pointer :: listNode(:)
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Field management
!
! Detect if node has equation defined in the field
!
! --------------------------------------------------------------------------------------------------
!
! In  field            : field
! In  nbNodeMesh       : number of nodes in mesh
! Ptr listNode         : pointer to list of nodes in mesh
!                       for iNode =  [1:nbNodeMesh]
!                           listNode[iNode] = 0 if node+component not present
!                           listNode[iNode] = 1 if node+component is present
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), pointer :: deeq(:) => null()
    character(len=19) :: nume_equa
    character(len=16) :: fieldSupp
    character(len=24) :: fieldRefe
    integer(kind=8) :: iEqua, iNodeMesh
    integer(kind=8) :: numeEqua, numeNode, numeCmp
    integer(kind=8) :: nbEquaChck, nbCmpOnNode, nbEqua
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM11_7')
    end if
!
! - Get parameters
!
    fieldRefe = field%fieldRefe
    nbEqua = field%nbEqua
    fieldSupp = field%fieldSupp
    ASSERT(fieldSupp .eq. 'NOEU')
!
! - Properties of field
!
    call dismoi('NUME_EQUA', fieldRefe, 'CHAM_NO', repk=nume_equa)
    call jeveuo(nume_equa(1:19)//'.DEEQ', 'L', vi=deeq)
    call jelira(fieldRefe(1:19)//'.VALE', 'LONMAX', nbEquaChck)
    ASSERT(nbEquaChck .eq. nbEqua)
!
! - Convert
!
    do iEqua = 1, nbEqua
! ----- Get current equation
        numeEqua = iEqua

! ----- Detect {Node,C mp} on current equation
        numeNode = deeq(2*(numeEqua-1)+1)
        numeCmp = deeq(2*(numeEqua-1)+2)

! ----- Physical node
        if (numeNode .gt. 0 .and. numeCmp .gt. 0) then
            listNode(numeNode) = listNode(numeNode)+1
        end if

! ----- Non-Physical node (Lagrange)
        if (numeNode .gt. 0 .and. numeCmp .lt. 0) then
            ASSERT(ASTER_FALSE)
        end if

! ----- Non-Physical node (Lagrange) - LIAISON_DDL
        if (numeNode .eq. 0 .and. numeCmp .eq. 0) then
            ASSERT(ASTER_FALSE)
        end if
    end do
!
! - Set value
!
    do iNodeMesh = 1, nbNodeMesh
        nbCmpOnNode = listNode(iNodeMesh)
        if (nbCmpOnNode .ne. 0) then
            listNode(iNodeMesh) = iNodeMesh
        end if
    end do
!
end subroutine
