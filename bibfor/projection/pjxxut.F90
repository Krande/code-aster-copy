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
subroutine pjxxut(projDime, typeSelect, &
                  entity1, entity2, &
                  nbCellSelect1, listCellSelect1, &
                  nbNodeSelect2, listNodeSelect2, &
                  mesh1, mesh2, &
                  nbCellType, cellListNume, cellListCode)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/pjnout.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
!
    character(len=2), intent(in) :: projDime
    character(len=*), intent(in) :: typeSelect
    character(len=8), intent(in) :: entity1, entity2
    integer(kind=8), intent(in) :: nbCellSelect1, listCellSelect1(*)
    integer(kind=8), intent(in) :: nbNodeSelect2, listNodeSelect2(*)
    character(len=8), intent(out) :: mesh1, mesh2
    integer(kind=8), intent(out) :: nbCellType, cellListNume(MT_NTYMAX)
    character(len=8), intent(out) :: cellListCode(MT_NTYMAX)
!
! --------------------------------------------------------------------------------------------------
!
! Projection of cells between two spaces
!
!   PREPARER LA LISTE DES MAILLES ET LES LISTES DE NOEUDS
!   UTILES A LA PROJECTION:
!
!   CETTE ROUTINE PRODUIT LES OBJETS SUIVANTS :
!    '&&PJXXCO.LIMA1' : NUMEROS DES MAILLES UTILES DE MOA1
!    '&&PJXXCO.LINO1' : NUMEROS DES NOEUDS UTILES DE MOA1
!    '&&PJXXCO.LINO2' : NUMEROS DES NOEUDS UTILES DE MOA2
!
!   M1 EST LE NOM DU MAILLAGE (OU DU MODELE) INITIAL
!   M2 EST LE NOM DU MAILLAGE (OU DU MODELE) FINAL
!
!   LES MAILLES UTILES DE MOA1 SONT CELLES QUI :
!      - SONT D'UN TYPE COHERENT AVEC DIM :
!             PAR EXEMPLE : '2D' -> TRIA/QUAD
!      - SONT PORTEUSES D'ELEMENTS FINIS (SI M1 EST UN MODELE)
!      - SONT INCLUSES DANS LIMA1 (SI MOCLE='PARTIE')
!
!   LES NOEUDS UTILES DE MOA1 SONT CEUX QUI SONT PORTES PAR LES
!   MAILLES UTILES DE MOA1
!
!   LES NOEUD UTILES DE MOA2 SONT CEUX QUI :
!      - SONT PORTES PAR LES MAILLES SUPPORTANT LES ELEMENTS FINIS
!        (SI M1 EST UN MODELE)
!      - SONT INCLUS DANS LINO2 (SI MOCLE='PARTIE')
!
!  SI MOCLE='TOUT' :
!     - ON NE SE SERT PAS DE NBMA1,LIMA1,NBNO2,LINO2
!
! --------------------------------------------------------------------------------------------------
!
! In  projDime         : dimension of projection 0D 1D 2D 3D
! In  typeSelect       : type of selection (all or restricted) TOUT or PARTIE
! In  entity1          : name of first entity (model or mesh)
! In  entity2          : name of second entity (model or mesh)
! In  nbCellSelect1    : number of cells to select in first entity if typeSelect = 'PARTIE'
! In  listCellSelect1  : list of cells to select in first entity if typeSelect = 'PARTIE'
! In  nbNodeSelect2    : number of nodes to select in second entity if typeSelect = 'PARTIE'
! In  listNodeSelect2  : list of nodes to select in second entity if typeSelect = 'PARTIE'
! Out mesh1            : name of mesh for first entity
! Out mesh2            : name of mesh for second entity
! Out nbCellType       : number of type of cells for this type of projDime
! Out cellListNume     : name of type of cells for this type of projDime
! Out cellListCode     : code of type of cells for this type of projDime
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model1, model2
    character(len=8) :: cellTypeName(MT_NTYMAX)
    integer(kind=8) :: nbCellModel, nbNode, nbNode1, nbNode2, nbCell1, nbCell2, nbNodeCount
    integer(kind=8) :: iCell1, iCellType, iNode, iCellModel, iCellSelect1, iNode2, iNodeSelect2
    integer(kind=8) :: nodeNume
    integer(kind=8) :: iad, ilcnx1, iexi
    character(len=19) :: ligrel
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: listCell1(:) => null(), listNode1(:) => null()
    integer(kind=8), pointer :: listNode2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    mesh1 = ' '
    mesh2 = ' '
    nbCellType = 0
    cellListNume = 0
    cellListCode = ' '

! - Detect model or mesh for first entity
    model1 = ' '
    call jeexin(entity1//'.MODELE    .NBNO', iexi)
    if (iexi .gt. 0) then
        model1 = entity1
        call dismoi('NOM_MAILLA', model1, 'MODELE', repk=mesh1)
    else
        model1 = ' '
        mesh1 = entity1
    end if

! - Detect model or mesh for second entity
    model2 = ' '
    call jeexin(entity2//'.MODELE    .NBNO', iexi)
    if (iexi .gt. 0) then
        model2 = entity2
        call dismoi('NOM_MAILLA', model2, 'MODELE', repk=mesh2)
        call pjnout(model2)
    else
        model2 = ' '
        mesh2 = entity2
    end if

! - Parameters about meshes
    call dismoi('NB_NO_MAILLA', mesh1, 'MAILLAGE', repi=nbNode1)
    call dismoi('NB_NO_MAILLA', mesh2, 'MAILLAGE', repi=nbNode2)
    call dismoi('NB_MA_MAILLA', mesh1, 'MAILLAGE', repi=nbCell1)
    call dismoi('NB_MA_MAILLA', mesh2, 'MAILLAGE', repi=nbCell2)

! - List of type of cells with this dimension
    if (projDime .eq. '0D') then
        nbCellType = 1
        cellTypeName(1) = 'POI1'
        cellListCode(1) = 'PO1'
    else if (projDime .eq. '1D') then
        nbCellType = 3
        cellTypeName(1) = 'SEG2'
        cellTypeName(2) = 'SEG3'
        cellTypeName(3) = 'SEG4'
        cellListCode(1) = 'SE2'
        cellListCode(2) = 'SE3'
        cellListCode(3) = 'SE4'

    else if (projDime .eq. '2D') then
        nbCellType = 6
        cellTypeName(1) = 'TRIA3'
        cellTypeName(2) = 'TRIA6'
        cellTypeName(3) = 'TRIA7'
        cellTypeName(4) = 'QUAD4'
        cellTypeName(5) = 'QUAD8'
        cellTypeName(6) = 'QUAD9'
        cellListCode(1) = 'TR3'
        cellListCode(2) = 'TR6'
        cellListCode(3) = 'TR7'
        cellListCode(4) = 'QU4'
        cellListCode(5) = 'QU8'
        cellListCode(6) = 'QU9'

    else if (projDime .eq. '3D') then
        nbCellType = 11
        cellTypeName(1) = 'TETRA4'
        cellTypeName(2) = 'TETRA10'
        cellTypeName(3) = 'PENTA6'
        cellTypeName(4) = 'PENTA15'
        cellTypeName(5) = 'PENTA18'
        cellTypeName(6) = 'HEXA8'
        cellTypeName(7) = 'HEXA20'
        cellTypeName(8) = 'HEXA27'
        cellTypeName(9) = 'PYRAM5'
        cellTypeName(10) = 'PYRAM13'
        cellTypeName(11) = 'HEXA9'
        cellListCode(1) = 'TE4'
        cellListCode(2) = 'T10'
        cellListCode(3) = 'PE6'
        cellListCode(4) = 'P15'
        cellListCode(5) = 'P18'
        cellListCode(6) = 'HE8'
        cellListCode(7) = 'H20'
        cellListCode(8) = 'H27'
        cellListCode(9) = 'PY5'
        cellListCode(10) = 'P13'
        cellListCode(11) = 'HE9'

    else
        ASSERT(ASTER_FALSE)
    end if

    do iCellType = 1, nbCellType
        call jenonu(jexnom('&CATA.TM.NOMTM', cellTypeName(iCellType)), cellListNume(iCellType))
    end do

! - Get list of cells for first entity
    call wkvect('&&PJXXCO.LIMA1', 'V V I', nbCell1, vi=listCell1)
    if (model1 .eq. ' ') then
        do iCell1 = 1, nbCell1
            listCell1(iCell1) = 1
        end do
    else
        call dismoi('NOM_LIGREL', model1, 'MODELE', repk=ligrel)
        call jeveuo(ligrel//'.TYFE', 'L', iad)
        call jelira(ligrel//'.TYFE', 'LONMAX', nbCellModel)
        do iCellModel = 1, nbCellModel
            if (zi(iad-1+iCellModel) .ne. 0) then
                listCell1(iCellModel) = 1
            end if
        end do
    end if

! - Select from type for first entity
    call jeveuo(mesh1//'.TYPMAIL', 'L', iad)
    do iCellType = 1, nbCellType
        do iCell1 = 1, nbCell1
            if (zi(iad-1+iCell1) .eq. cellListNume(iCellType)) then
                listCell1(iCell1) = listCell1(iCell1)+1
            end if
        end do
    end do

! - Final list of 'good' cells for first entity
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .eq. 1) then
            listCell1(iCell1) = 0
        else if (listCell1(iCell1) .eq. 2) then
            listCell1(iCell1) = 1
        else if (listCell1(iCell1) .gt. 2) then
            ASSERT(ASTER_FALSE)
        end if
    end do

! - Selection of cells for first entity when list is given
    if (typeSelect .eq. 'PARTIE') then
        do iCellSelect1 = 1, nbCellSelect1
            listCell1(listCellSelect1(iCellSelect1)) = 2*listCell1(listCellSelect1(iCellSelect1))
        end do
        do iCell1 = 1, nbCell1
            listCell1(iCell1) = listCell1(iCell1)/2
        end do
    end if

! - Get list of 'good' nodes for first entity
    call wkvect('&&PJXXCO.LINO1', 'V V I', nbNode1, vi=listNode1)
    call jeveuo(mesh1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            nbNode = zi(ilcnx1+iCell1)-zi(ilcnx1-1+iCell1)
            do iNode = 1, nbNode
                nodeNume = connex(1+zi(ilcnx1-1+iCell1)-2+iNode)
                listNode1(nodeNume) = 1
            end do
        end if
    end do

! - Get list of 'good' nodes for second entity
    call wkvect('&&PJXXCO.LINO2', 'V V I', nbNode2, vi=listNode2)
    if (model2 .ne. ' ') then
        call jeveuo(model2//'.NOEUD_UTIL', 'L', iad)
        if (typeSelect .eq. 'TOUT') then
            do iNode2 = 1, nbNode2
                if (zi(iad-1+iNode2) .ne. 0) then
                    listNode2(iNode2) = 1
                end if
            end do
        else if (typeSelect .eq. 'PARTIE') then
            do iNodeSelect2 = 1, nbNodeSelect2
                if (zi(iad-1+listNodeSelect2(iNodeSelect2)) .ne. 0) then
                    listNode2(listNodeSelect2(iNodeSelect2)) = 1
                end if
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        if (typeSelect .eq. 'TOUT') then
            do iNode2 = 1, nbNode2
                listNode2(iNode2) = 1
            end do
        else if (typeSelect .eq. 'PARTIE') then
            do iNodeSelect2 = 1, nbNodeSelect2
                listNode2(listNodeSelect2(iNodeSelect2)) = 1
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Stop if there node nodes in second entity
    nbNodeCount = 0
    do iNode2 = 1, nbNode2
        if (listNode2(iNode2) .gt. 0) then
            nbNodeCount = nbNodeCount+1
        end if
    end do
    if (nbNodeCount .eq. 0) then
        call utmess('F', 'PROJECTION4_54')
    end if
!
    call jedema()
end subroutine
