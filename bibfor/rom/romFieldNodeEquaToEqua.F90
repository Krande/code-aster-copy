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
subroutine romFieldNodeEquaToEqua(fieldA, fieldB, nbNodeMesh, listNode, equaAToB)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/romFieldNodeFromEqua.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_Field), intent(in) :: fieldA, fieldB
    integer(kind=8), intent(in) :: nbNodeMesh
    integer(kind=8), pointer :: listNode(:), equaAToB(:)
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Create map for equation numbering from domain A to domain B - Nodal field
!
! --------------------------------------------------------------------------------------------------
!
! The object listNode has been created before to be efficient.
!
! In  fieldA           : field (representative) on domain A
! In  fieldB           : field (representative) on domain B
! In  nbNodeMesh       : number of nodes in mesh
! Ptr listNode         : pointer to list of nodes in mesh for selection
!                          for iNode =  [1:nbNodeMesh]
!                              listNode[iNode] = 0 if node+component not present from fieldA
!                              listNode[iNode] = 1 if node+component is present from fieldA
! Ptr equaAToB         : pointer to map for equation numbering from domain A to domain B
!                          for iEquaA =  [1:nbEquaA]
!                              equaAToB[iEquaA] = 0 => this equation is not in domain B
!                              equaAToB[iEquaA] = numeEquaB => this equation is in domain B
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: iLigrMesh = 1
    integer(kind=8), parameter :: nbEcMax = 10
    character(len=24) :: fieldRefeA, fieldRefeB
    character(len=16) :: fieldSuppA, fieldSuppB
    integer(kind=8) :: nb_ec, nbEquaA, nbEquaB
    integer(kind=8) :: iEquaA, physNumeA, physNumeB, numeDofB
    integer(kind=8) :: numeNode, numeCmp, numeEquaA, numeEquaB
    character(len=19) :: profChnoA, profChnoB
    integer(kind=8), pointer :: deeqA(:) => null(), deeqB(:) => null()
    integer(kind=8), pointer :: nueqB(:) => null()
    integer(kind=8), pointer :: prnoB(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM11_37')
    end if
!
! - Parameters on fields
!
    fieldRefeA = fieldA%fieldRefe
    fieldSuppA = fieldA%fieldSupp
    nbEquaA = fieldA%nbEqua
    ASSERT(fieldSuppA .eq. 'NOEU')
    fieldRefeB = fieldB%fieldRefe
    fieldSuppB = fieldB%fieldSupp
    nbEquaB = fieldB%nbEqua
    ASSERT(fieldSuppB .eq. 'NOEU')
!
! - Access to numbering for domain A
!
    call dismoi('NUME_EQUA', fieldRefeA, 'CHAM_NO', repk=profChnoA)
    call jeveuo(profChnoA(1:19)//'.DEEQ', 'L', vi=deeqA)
!
! - Access to numbering for domain B
!
    call dismoi('NUME_EQUA', fieldRefeB, 'CHAM_NO', repk=profChnoB)
    call jeveuo(profChnoB(1:19)//'.DEEQ', 'L', vi=deeqB)
    call jeveuo(jexnum(profChnoB(1:19)//'.PRNO', iLigrMesh), 'L', vi=prnoB)
    call jeveuo(profChnoB(1:19)//'.NUEQ', 'L', vi=nueqB)
!
! - Get informations about physical quantity
!
    physNumeA = 0
    nb_ec = 0
    call dismoi('NUM_GD', fieldRefeA, 'CHAM_NO', repi=physNumeA)
    ASSERT(physNumeA .ne. 0)
    call dismoi('NUM_GD', fieldRefeB, 'CHAM_NO', repi=physNumeB)
    ASSERT(physNumeA .eq. physNumeB)
    nb_ec = nbec(physNumeA)
    ASSERT(nb_ec .le. nbEcMax)
!
! - Detect nodes on all mesh with equation from reduced domain
!
    call romFieldNodeFromEqua(fieldB, nbNodeMesh, listNode)
!
! - Select equations
!
    do iEquaA = 1, nbEquaA
        numeEquaA = iEquaA
! ----- Get equation information in domain A
        numeNode = deeqA(2*(numeEquaA-1)+1)
        numeCmp = deeqA(2*(numeEquaA-1)+2)
! ----- For physical node
        if (numeNode .gt. 0 .and. numeCmp .gt. 0) then
            if (listNode(numeNode) .eq. 0) then
! ------------- This node is NOT affected in domain B
                equaAToB(numeEquaA) = 0
            else
! ------------- This node is in affected in domain B
                numeDofB = prnoB((nb_ec+2)*(numeNode-1)+1)-1
! ------------- Index of equation in domain B
                numeEquaB = nueqB(numeDofB+numeCmp)
                ASSERT(deeqB(2*(numeEquaB-1)+1) .eq. numeNode)
                ASSERT(deeqB(2*(numeEquaB-1)+2) .eq. numeCmp)
! ------------- Save index of equation in domain B
                equaAToB(numeEquaA) = numeEquaB
            end if
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
! - Reinitialization of list of nodes
!
    listNode = 0
!
end subroutine
