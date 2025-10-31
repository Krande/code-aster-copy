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
subroutine select_dof(listEqua_, tablEqua_, tablCmp_, &
                      numeDofZ_, fieldNodeZ_, &
                      nbNodeToSelect_, listNodeToSelect_, &
                      nbCmpToSelect_, listCmpToSelect_)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/select_dof_gene.h"
!
    integer(kind=8), pointer, optional :: listEqua_(:), tablEqua_(:, :)
    integer(kind=8), pointer, optional :: tablCmp_(:)
    character(len=*), optional, intent(in) :: numeDofZ_, fieldNodeZ_
    integer(kind=8), optional, intent(in) :: nbNodeToSelect_
    integer(kind=8), pointer, optional :: listNodeToSelect_(:)
    integer(kind=8), optional, intent(in) :: nbCmpToSelect_
    character(len=8), pointer, optional :: listCmpToSelect_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Utility for nodal field
!
! Select dof from list of nodes and components - Only on mesh
!
! --------------------------------------------------------------------------------------------------
!
! Select output: listEqua/tablEqua/tablCmp
! Select input:
!    On numbering of nodal field (numeDof / fieldNode)
!    On which nodes (nbNodeToSelect / listNodeToSelect)
!                   if none => all nodes in mesh
!    On which components (nbCmpToSelect / listCmpToSelect)
!                   if none => all components
!
! Ptr listEqua             : pointer to all equations [1:nbEqua] with
!                       for iEqua =  [1:nbEqua]
!                           listEqua[iEqua] = 0 if node+component not present
!                           listEqua[iEqua] = 1 if node+component is present
! Ptr tablEqua             : pointer to all equations and all components
!                            [1:nbEqua, 1:nbCmpToSelect] with
!                       for iEqua = [1:nbEqua]
!                           for iCmp = [1:nbCmpToSelect]
!                               tablEqua[iEqua,iCmp] = 0 if node+component not present
!                               tablEqua[iEqua,iCmp] = 1 if node+component is present
! Ptr tablCmp_             : pointer to all components [1:nbCmpToSelect]
!                       for iCmp = [1:nbCmpToSelect]
!                           tablCmp_[iCmp] = 0   if node+component not present
!                           tablCmp_[iCmp] = iEqua if node+component present
! In  numeDof              : name of numbering (NUME_DDL)
! In  fieldNode            : name of nodal field (CHAMNO)
! In  nbNodeToSelect       : number of nodes to select (if none => all nodes in mesh)
! Ptr listNodeToSelect     : pointer to the list of nodes (absolute index in mesh)
! In  nbCmpToSelect        : number of components (if none => all components)
! Ptr listCmpToSelect      : pointer to the list of components (name)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: iLigrMesh = 1
    integer(kind=8), parameter :: nbEcMax = 11
    integer(kind=8) :: physDesc(nbEcMax)
    character(len=24) :: liliName
    character(len=8) :: cmpName, mesh
    character(len=19) :: nume_equa, nume_equa_gene, fieldNode, numeEqul
    character(len=14) :: numeDof
    integer(kind=8) :: iexi
    aster_logical :: lnume_equa_gene
    aster_logical :: lMatrDist
    integer(kind=8) :: nodeNume, physNume, prnoLength
    integer(kind=8) :: nbNodeToSelect
    integer(kind=8) :: numeEqua, iNode, iCmp, dofNume, numeEquaL, numeCmp, iEc, nbNodeMesh
    integer(kind=8) :: nb_ec, nbCmpToSelect, physNbCmp, nbCmpNode
    integer(kind=8), pointer :: cmpSelect(:) => null()
    integer(kind=8), pointer :: nodeSelect(:) => null()
    integer(kind=8), pointer :: prno(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: nugl(:) => null()
    character(len=8), pointer :: physCataName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    numeDof = ' '
    fieldNode = ' '
    lMatrDist = ASTER_FALSE
!
! - Check output parameters
!
    ASSERT(EXCLUS2(tablEqua_, tablCmp_))
    ASSERT(EXCLUS2(tablEqua_, listEqua_))
    ASSERT(EXCLUS2(listEqua_, tablCmp_))
!
! - Check input parameters
!
    if (present(tablCmp_)) then
        ASSERT(present(nbNodeToSelect_))
        ASSERT(present(listNodeToSelect_))
        ASSERT(nbNodeToSelect_ .eq. 1)
    end if
!
! - Get NUME_EQUA
!
    nume_equa = ' '
    if (present(numeDofZ_)) then
        numeDof = numeDofZ_
        ASSERT(.not. present(fieldNodeZ_))
        call dismoi('NUME_EQUA', numeDof, 'NUME_DDL', repk=nume_equa)
    elseif (present(fieldNodeZ_)) then
        fieldNode = fieldNodeZ_
        ASSERT(.not. present(numeDofZ_))
        call dismoi('NUME_EQUA', fieldNode, 'CHAM_NO', repk=nume_equa)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Check if nume_ddl is correct (Distributed matrix)
!
    if (present(numeDofZ_)) then
        numeEqul = numeDof//'.NUML'
        call jeexin(numeEqul(1:19)//'.NUGL', iexi)
        lMatrDist = iexi .ne. 0
        if (lMatrDist) then
            call jeveuo(numeEqul(1:19)//'.NUGL', 'L', vi=nugl)
        end if
    end if
!
! - Get informations about physical quantity
!
    physNume = 0
    nb_ec = 0
    if (present(numeDofZ_)) then
        call dismoi('NUM_GD_SI', numeDof, 'NUME_DDL', repi=physNume)
    elseif (present(fieldNodeZ_)) then
        call dismoi('NUM_GD', fieldNode, 'CHAM_NO', repi=physNume)
    else
        ASSERT(ASTER_FALSE)
    end if
    ASSERT(physNume .ne. 0)
    nb_ec = nbec(physNume)
    ASSERT(nb_ec .le. nbEcMax)
!
! - Access to catalog
!
    call jelira(jexnum('&CATA.GD.NOMCMP', physNume), 'LONMAX', physNbCmp)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', physNume), 'L', vk8=physCataName)
!
! - Select number of components
!
    if (present(listCmpToSelect_)) then
        ASSERT(present(nbCmpToSelect_))
        nbCmpToSelect = nbCmpToSelect_
    else
        nbCmpToSelect = physNbCmp
    end if
!
! - Select components
!
    AS_ALLOCATE(vi=cmpSelect, size=physNbCmp)
    do iCmp = 1, nbCmpToSelect
        if (present(listCmpToSelect_)) then
            cmpName = listCmpToSelect_(iCmp)
        else
            cmpName = physCataName(iCmp)
        end if
        numeCmp = indik8(physCataName, cmpName, 1, physNbCmp)
        if (numeCmp .ne. 0) then
            cmpSelect(numeCmp) = iCmp
        end if
    end do
!
! - NUME_EQUA or NUME_EQUA_GENE ?
!
    call jeexin(nume_equa//'.DESC', iexi)
    lnume_equa_gene = (iexi .gt. 0)
    if (lnume_equa_gene) then
        nume_equa_gene = nume_equa
        call select_dof_gene(nume_equa_gene, nbCmpToSelect, physCataName, &
                             listCmpToSelect_, listEqua_, tablEqua_)
        goto 99
    end if
!
! - Get mesh
!
    mesh = ' '
    if (present(numeDofZ_)) then
        call dismoi('NOM_MAILLA', numeDof, 'NUME_DDL', repk=mesh)
    elseif (present(fieldNodeZ_)) then
        call dismoi('NOM_MAILLA', fieldNode, 'CHAM_NO', repk=mesh)
    else
        ASSERT(ASTER_FALSE)
    end if
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbNodeMesh)
!
! - Select number of nodes
!
    nbNodeToSelect = 0
    if (present(listNodeToSelect_)) then
        nbNodeToSelect = nbNodeToSelect_
    else
        nbNodeToSelect = nbNodeMesh
    end if
!
! - Select nodes
!
    AS_ALLOCATE(vi=nodeSelect, size=nbNodeToSelect)
    if (present(listNodeToSelect_)) then
        do iNode = 1, nbNodeToSelect
            nodeNume = listNodeToSelect_(iNode)
            nodeSelect(iNode) = nodeNume
        end do
    else
        do iNode = 1, nbNodeToSelect
            nodeNume = iNode
            nodeSelect(iNode) = nodeNume
        end do
    end if
!
! - Some checks
!
    call jenuno(jexnum(nume_equa(1:19)//'.LILI', iLigrMesh), liliName)
    ASSERT(liliName .eq. '&MAILLA')
    call jelira(jexnum(nume_equa(1:19)//'.PRNO', iLigrMesh), 'LONMAX', prnoLength)
    ASSERT(prnoLength/(nb_ec+2) .eq. nbNodeMesh)
!
! - Get objects
!
    call jeveuo(jexnum(nume_equa(1:19)//'.PRNO', iLigrMesh), 'L', vi=prno)
    call jeveuo(nume_equa(1:19)//'.NUEQ', 'L', vi=nueq)
!
! - Loop on nodes
!
    if (nbNodeToSelect .ne. 0) then
        do iNode = 1, nbNodeToSelect
            nodeNume = nodeSelect(iNode)

! --------- Parameters of current node
            dofNume = prno((nb_ec+2)*(nodeNume-1)+1)-1
            nbCmpNode = prno((nb_ec+2)*(nodeNume-1)+2)

! --------- Vector containing active components on current node
            physDesc = 0
            if (nbCmpNode .ne. 0) then
                do iEc = 1, nb_ec
                    physDesc(iEc) = prno((nb_ec+2)*(nodeNume-1)+2+iEc)
                end do
            end if

! --------- Loop on components to seek
            do iCmp = 1, physNbCmp
                if (exisdg(physDesc, iCmp)) then
                    dofNume = dofNume+1
                    numeCmp = cmpSelect(iCmp)
                    if (numeCmp .ne. 0) then
                        numeEqua = nueq(dofNume)
                        if (lMatrDist) then
                            numeEquaL = nugl(numeEqua)
                        else
                            numeEquaL = numeEqua
                        end if
                        if (present(tablCmp_)) then
                            tablCmp_(numeCmp) = numeEquaL
                        elseif (present(listEqua_)) then
                            listEqua_(numeEquaL) = 1
                        elseif (present(tablEqua_)) then
                            tablEqua_(numeEquaL, numeCmp) = 1
                        end if
                    end if
                end if
            end do
        end do
    end if
!
99  continue
!
! - Clean
!
    AS_DEALLOCATE(vi=cmpSelect)
    AS_DEALLOCATE(vi=nodeSelect)
!
end subroutine
