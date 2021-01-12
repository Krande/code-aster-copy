! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1501
!
subroutine op0167()
!
use mesh_module, only : checkInclude, createNameOfCell,&
                        getCellOptionForName, getNodeOptionForName
use SolidShell_Mesh_module, only : orieHexa9
use crea_maillage_module
!
implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/chckma.h"
#include "asterfort/chcoma.h"
#include "asterfort/chcomb.h"
#include "asterfort/cm_dclac.h"
#include "asterfort/cm1518.h"
#include "asterfort/cm2027.h"
#include "asterfort/cmhho.h"
#include "asterfort/cmcovo.h"
#include "asterfort/cmcrea.h"
#include "asterfort/cmlqlq.h"
#include "asterfort/cmmoma.h"
#include "asterfort/cmqlql.h"
#include "asterfort/cmqutr.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/eclpgm.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infoma.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rdtmai.h"
#include "asterfort/getelem.h"
#include "asterfort/reliem.h"
#include "asterfort/getnode.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! CREA_MAILLAGE
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nori, ntab, n1
    integer :: k, iret, iqtr
    real(kind=8) :: epais, shrink, lonmin
    character(len=4) :: answer
    character(len=4) :: cdim
    character(len=8) :: meshIn, meshOut, model, geofi
    character(len=8) :: nomori, plan, trans, knume
    character(len=16) :: option, keywfact
    character(len=16) :: kbi1, kbi2
    character(len=19) :: table, ligrel
    character(len=19), parameter :: k19void = ' '
    integer :: nbMeshIn
    character(len=24), parameter :: crgrnu = '&&OP0167.CR_GR.NUM'
    character(len=24), parameter :: crgrno = '&&OP0167.CR_GR.NOM'
    character(len=24) :: grCellName, grNodeName
    integer :: jvConnexIn, jvConnexOut, jvGeofi
    integer :: connexLength, cellShift, listCreaLength
    integer :: cellNumeIn, cellNumeOut, cellTypeIn, cellTypeToModify
    integer :: nodeNumeIn, nodeNumeOut
    integer :: iCell, iNode, iCellModi, iGrCell, iGrNode
    integer :: iocc, nbOcc
    character(len=24), parameter :: jvCellNume = '&&OP0167.LISTCELL'
    character(len=24), parameter :: jvNodeNume = '&&OP0167.LISTNODE'
    integer :: nbCellIn, nbCellOut, nbCell, nbCellCrea, nbCellModi, nbCellType, nbCellAddPoi1
    integer :: nbNodeInCellOut, nbNodeInCellIn
    integer :: nbGrCellFromCreaCell, nbGrCellFromCreaPoi1, nbGrCellIn, nbGrCellOut
    integer :: nbGrNodeIn, nbGrNodeOut
    integer :: nbCellInGrOut, nbCellInGrIn, nbNodeInGrOut, nbNodeInGrIn
    integer :: nbNodeCrea, nbNodeIn, nbNodeOut, nbNode
    integer :: nbField
    integer :: nbOccDecoupeLac, nbOccEclaPg, nbGeomFibre, nbOccCreaFiss, nbOccLineQuad
    integer :: nbOccQuadLine, nbOccModiMaille, nbOccCoquVolu, nbOccRestreint, nbOccRepere
    integer :: iOccQuadTria, iad
    integer :: nbOccCreaPoi1, nbOccCreaMaille, nbOccModiHHO, nbOccCoqueSolide
    aster_logical :: lpb
    character(len=8) :: cellNameIn, cellNameOut, nodeNameIn, nodeNameOut
    aster_logical :: lPrefCellName, lPrefCellNume, lPrefNodeName, lPrefNodeNume
    integer :: prefCellNume, prefNodeNume, prefNume
    character(len=8) :: prefCellName, prefNodeName
    integer, pointer :: modiCellNume(:) => null(), modiCellType(:) => null()
    integer, pointer :: listCellNume(:) => null(), listNodeNume(:) => null()
    character(len=16), pointer :: listField(:) => null()
    integer, pointer :: connexAdr(:) => null()
    integer, pointer :: nbNodeByCellOut(:) => null()
    integer, pointer :: allCellNume(:) => null()
    character(len=24), pointer :: creaCellOccGrName(:) => null()
    character(len=8), pointer :: creaCellName(:) => null()
    integer, pointer :: creaCellNume(:) => null()
    integer, pointer :: creaCellOccNb(:) => null()
    character(len=8), pointer :: creaNodeName(:) => null()
    integer, pointer :: creaNodeNume(:) => null()
    character(len=8), pointer :: creaPoi1Name(:) => null(), creaGrPoi1CellName(:) => null()
    aster_logical, pointer :: creaPoi1Flag(:) => null()
    integer, pointer :: creaGrPoi1NbCell(:) => null()
    character(len=24), pointer :: creaGrPoi1GrName(:) => null()
    integer, pointer :: cellInGrIn(:) => null(), cellInGrOut(:) => null()
    integer, pointer :: nodeInGrIn(:) => null(), nodeInGrOut(:) => null()
    integer, pointer :: meshDimeIn(:) => null(), meshDimeOut(:) => null()
    integer, pointer :: meshTypmailIn(:) => null(), meshTypmailOut(:) => null()
    character(len=24), pointer :: meshRefeOut(:) => null()
    real(kind=8), pointer :: meshValeIn(:) => null(), meshValeOut(:) => null()
    type(Mmesh) :: meshSolidShell
    character(len=8) :: convType(2)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
! - Initializations
!
    nbNodeCrea           = 0
    nbCellCrea           = 0
    nbCellAddPoi1        = 0
    nbGrCellFromCreaCell = 0
    nbGrCellFromCreaPoi1 = 0
    nbGrCellIn           = 0
    nbGrNodeIn           = 0
!
! - Check Includes
!
    call checkInclude()
!
! - Keywords
!
    call getfac('ECLA_PG', nbOccEclaPg)
    call getvid(' ', 'GEOM_FIBRE', scal=geofi, nbret=nbGeomFibre)
    call getfac('CREA_FISS', nbOccCreaFiss)
    call getfac('LINE_QUAD', nbOccLineQuad)
    call getfac('QUAD_LINE', nbOccQuadLine)
    call getfac('MODI_MAILLE', nbOccModiMaille)
    call getfac('COQU_VOLU', nbOccCoquVolu)
    call getfac('RESTREINT', nbOccRestreint)
    call getfac('REPERE', nbOccRepere)
    call getfac('DECOUPE_LAC', nbOccDecoupeLac)
    call getfac('CREA_MAILLE', nbOccCreaMaille)
    call getfac('CREA_POI1', nbOccCreaPoi1)
    call getfac('MODI_HHO', nbOccModiHHO)
    call getfac('COQUE_SOLIDE', nbOccCoqueSolide)
!
! - Main datastructure
!
    call getres(meshOut, kbi1, kbi2)
    call getvid(' ', 'MAILLAGE', scal = meshIn, nbret = nbMeshIn)
    if (nbMeshIn .ne. 0) then
        if (isParallelMesh(meshIn)) then
            call utmess('F', 'MESH1_22')
        end if
        call jeveuo(meshIn//'.DIME', 'L', vi = meshDimeIn)
        nbNodeIn   = meshDimeIn(1)
        nbCellIn   = meshDimeIn(3)
        call jeexin(meshIn//'.GROUPEMA', iret)
        if (iret .ne. 0) then
            call jelira(meshIn//'.GROUPEMA', 'NOMUTI', nbGrCellIn)
        endif
        call jeexin(meshIn//'.GROUPENO', iret)
        if (iret .ne. 0) then
            call jelira(meshIn//'.GROUPENO', 'NOMUTI', nbGrNodeIn)
        endif
        call jeveuo(meshIn//'.TYPMAIL', 'L', vi = meshTypmailIn)
    endif
!
! - MAILLAGE is required except for GEOM_FIBRE and ECLA_PG
!
    if (nbMeshIn .eq. 0 .and. nbGeomFibre .eq. 0 .and. nbOccEclaPg .eq. 0) then
        call utmess('F', 'MESH1_15')
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "ECLA_PG"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccEclaPg .gt. 0) then
        keywfact = 'ECLA_PG'
        ASSERT(nbOccEclaPg .eq. 1)
        call getvid(keywfact, 'MODELE', iocc=1, scal=model)
        call getvr8(keywfact, 'SHRINK', iocc=1, scal=shrink)
        call getvr8(keywfact, 'TAILLE_MIN', iocc=1, scal=lonmin)
        call getvtx(keywfact, 'NOM_CHAM', iocc=1, nbval=0, nbret=nbField)
        if (nbField .lt. 0) then
            nbField = -nbField
            AS_ALLOCATE(vk16 = listField, size = nbField)
            call getvtx(keywfact, 'NOM_CHAM', iocc=1, nbval=nbField, vect=listField)
        else
            AS_ALLOCATE(vk16 = listField, size = 1)
        endif
        call exlima(keywfact, 1, 'V', model, ligrel)
        call eclpgm(meshOut, model  , k19void  , ligrel, shrink,&
                    lonmin , nbField, listField)
        AS_DEALLOCATE(vk16 = listField)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "GEOM_FIBRE"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbGeomFibre .ne. 0) then
        call jeveuo(geofi//'.GFMA', 'L', jvGeofi)
        call copisd('MAILLAGE', 'G', zk8(jvGeofi), meshOut)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "CREA_FISS"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccCreaFiss .ne. 0) then
        call cmcrea(meshIn, meshOut, nbOccCreaFiss)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "LINE_QUAD"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccLineQuad .gt. 0) then
        ASSERT(nbOccLineQuad .eq. 1)
        call jeexin(meshIn//'.NOMACR', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_7')
        endif
        call jeexin(meshIn//'.ABSC_CURV', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_8')
        endif
        keywfact = 'LINE_QUAD'
        call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNodeName)
        call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNodeNume)
        call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
        if (nbCell .ne. nbCellIn) then
            call utmess('A', 'MESH1_4', sk=keywfact)
        endif
        call jeveuo(jvCellNume, 'L', vi = listCellNume)
        call cmlqlq(meshIn, meshOut, nbCell, listCellNume, prefNodeName, prefNodeNume)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "PENTA15_18","HEXA20_27"
!
! --------------------------------------------------------------------------------------------------
!
    do k = 1, 2
        if (k .eq. 1) keywfact='HEXA20_27'
        if (k .eq. 2) keywfact='PENTA15_18'
        call getfac(keywfact, nbOcc)
        if (nbOcc .gt. 0) then
!
            lpb=.false.
            if (keywfact .eq. 'HEXA20_27') then
                call dismoi('EXI_PENTA15', meshIn, 'MAILLAGE', repk=answer)
                if (answer .eq. 'OUI') lpb=.true.
                call dismoi('EXI_PYRAM13', meshIn, 'MAILLAGE', repk=answer)
                if (answer .eq. 'OUI') lpb=.true.
            else if (keywfact.eq.'PENTA15_18') then
                call dismoi('EXI_HEXA20', meshIn, 'MAILLAGE', repk=answer)
                if (answer .eq. 'OUI') lpb=.true.
                call dismoi('EXI_PYRAM13', meshIn, 'MAILLAGE', repk=answer)
                if (answer .eq. 'OUI') lpb=.true.
            endif
            if (lpb) then
                call utmess('A', 'MESH1_11', sk=keywfact)
            endif
!
            call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNodeName)
            call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNodeNume)
!
            call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
            call jeveuo(jvCellNume, 'L', vi = listCellNume)
            if (nbCell .ne. nbCellIn) then
                call utmess('A', 'MESH1_4', sk=keywfact)
            endif
!
            if (keywfact .eq. 'HEXA20_27') then
                call cm2027(meshIn, meshOut, nbCell, listCellNume, prefNodeName, prefNodeNume)
            else if (keywfact.eq.'PENTA15_18') then
                call cm1518(meshIn, meshOut, nbCell, listCellNume, prefNodeName, prefNodeNume)
            endif
            goto 350
        endif
    end do
!
! --------------------------------------------------------------------------------------------------
!
!   For "COQUE_SOLIDE"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccCoqueSolide .gt. 0) then
        call jeexin(meshIn//'.NOMACR', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_7')
        endif
        call jeexin(meshIn//'.ABSC_CURV', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_8')
        endif
        keywfact = 'COQUE_SOLIDE'

! ----- Create mesh to convert
        call meshSolidShell%init(meshIn)

! ----- Add conversions
        convType = ["HEXA8", "HEXA9"]
        call meshSolidShell%converter%add_conversion(convType(1), convType(2))
        convType = ["PENTA6", "PENTA7"]
        call meshSolidShell%converter%add_conversion(convType(1), convType(2))

        do iocc = 1, nbOccCoqueSolide

! --------- Get parameters
            call getelem(meshIn, keywfact, iocc, 'F', jvCellNume, nbCell)
            call jeveuo(jvCellNume, 'L', vi = listCellNume)
            call getvtx(keywfact, 'PREF_NOEUD', iocc = iocc, scal = prefNodeName)
            call getvis(keywfact, 'PREF_NUME' , iocc = iocc, scal = prefNodeNume)

! --------- Convert cells
            call meshSolidShell%convert_cells(nbCell, listCellNume, prefNodeName, prefNodeNume)

! --------- Orient HEXA9
            call orieHexa9(iOcc, nbCell, listCellNume, meshIn)

            call jedetr(jvCellNume)
        end do
! ----- Copy mesh
        call meshSolidShell%copy_mesh(meshOut)
        call meshSolidShell%create_joints(meshOut)
        call meshSolidShell%clean()
        goto 350
    endif

!
! --------------------------------------------------------------------------------------------------
!
!   For "QUAD_LINE"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccQuadLine .gt. 0) then
        ASSERT(nbOccQuadLine .eq. 1)
        call jeexin(meshIn//'.NOMACR', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_7')
        endif
        call jeexin(meshIn//'.ABSC_CURV', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_8')
        endif
        keywfact = 'QUAD_LINE'
        call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
        if (nbCell .ne. nbCellIn) then
            call utmess('A', 'MESH1_4', sk=keywfact)
        endif
        call jeveuo(jvCellNume, 'L', vi = listCellNume)
        call cmqlql(meshIn, meshOut, nbCell, listCellNume)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "MODI_HHO"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccModiHHO .gt. 0) then
        ASSERT(nbOccModiHHO .eq. 1)
        call jeexin(meshIn//'.NOMACR', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_7')
        endif
        call jeexin(meshIn//'.ABSC_CURV', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_8')
        endif
        keywfact = 'MODI_HHO'
        call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
        if (nbCell .ne. nbCellIn) then
            call utmess('A', 'MESH1_4', sk=keywfact)
        endif
        call jeveuo(jvCellNume, 'L', vi = listCellNume)
        call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNodeName)
        call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNodeNume)
        call cmhho(meshIn, meshOut, nbCell, listCellNume, prefNodeName, prefNodeNume)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "MODI_MAILLE", OPTION "QUAD_TRIA3"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccModiMaille .gt. 0) then
        keywfact = 'MODI_MAILLE'
        iqtr = 0
        do iocc = 1, nbOccModiMaille
            call getvtx(keywfact, 'OPTION', iocc=iocc, scal=option)
            if (option .eq. 'QUAD_TRIA3') then
                iqtr         = iqtr + 1
                iOccQuadTria = iocc
            endif
        end do
        if (iqtr .gt. 1) then
            call utmess('F', 'MESH1_9')
        endif
        if (iqtr .eq. 1) then
            call dismoi('EXI_TRIA6', meshIn, 'MAILLAGE', repk=answer)
            if (answer .eq. 'OUI') then
                call utmess('A', 'MESH1_10')
            endif
            call getvtx(keywfact, 'PREF_MAILLE', iocc=iOccQuadTria, scal=prefCellName)
            call getvis(keywfact, 'PREF_NUME', iocc=iOccQuadTria, scal=prefCellNume)
            call getelem(meshIn, keywfact, iOccQuadTria, 'F', jvCellNume, nbCell)
            if (nbCell .ne. nbCellIn) then
                call utmess('A', 'MESH1_4', sk=keywfact)
            endif
            call jeveuo(jvCellNume, 'L', vi = listCellNume)
            call cmqutr('G', meshIn, meshOut, nbCell, listCellNume,&
                        prefCellName, prefCellNume)
            goto 350
        endif
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "COQU_VOLU"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccCoquVolu .ne. 0) then
        ASSERT(nbOccCoquVolu .eq. 1)
        keywfact = 'COQU_VOLU'
        call getvr8(keywfact, 'EPAIS', iocc=1, scal=epais)
        call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNodeName)
        call getvtx(keywfact, 'PREF_MAILLE', iocc=1, scal=prefCellName)
        call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNume)
        call getvtx(keywfact, 'PLAN', iocc=1, scal=plan)
        if (plan .eq. 'MOY') then
            trans = 'INF'
            call getvtx(keywfact, 'TRANSLATION', iocc=1, scal=trans)
        endif
        call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
        call jeveuo(jvCellNume, 'L', vi = listCellNume)
        call cmcovo(meshIn, meshOut, nbCell, listCellNume, prefNodeName,&
                    prefCellName, prefNume, epais, plan, trans)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "RESTREINT"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccRestreint .ne. 0) then
        call rdtmai(meshIn, meshOut, 'G', meshOut//'.CRNO', meshOut// '.CRMA', 'G', 0, [0])
        call chckma(meshOut, 1.0d-03)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "MODI_MAILLE"
!
! --------------------------------------------------------------------------------------------------
!
    nbNodeCrea = 0
    nbCellModi = 0
    if (nbOccModiMaille .ne. 0) then
        keywfact = 'MODI_MAILLE'
        AS_ALLOCATE(vi = modiCellNume, size = nbCellIn)
        AS_ALLOCATE(vi = modiCellType, size = nbCellIn)
        AS_ALLOCATE(vk8 = creaNodeName, size = nbCellIn)
        AS_ALLOCATE(vi  = creaNodeNume, size = nbCellIn)

        iad = 1
        do iocc = 1, nbOccModiMaille
! --------- How to transform cell ?
            call getvtx(keywfact, 'OPTION', iocc=iocc, scal=option)
            if (option .eq. 'TRIA6_7') then
                cellTypeToModify = MT_TRIA6
            else if (option .eq. 'QUAD8_9') then
                cellTypeToModify = MT_QUAD8
            else if (option .eq. 'SEG3_4') then
                cellTypeToModify = MT_SEG3
            else
                ASSERT(ASTER_FALSE)
            endif

! --------- Get options from user for name of new nodes
            call getNodeOptionForName(keywfact     , iocc,&
                                      lPrefNodeName, lPrefNodeNume,&
                                      prefNodeName , prefNodeNume)

! --------- Get list of cells to modify
            call getelem(meshIn, keywfact, iocc, 'F', jvCellNume, nbCell)
            call jeveuo(jvCellNume, 'L', vi = listCellNume)

! --------- Count number of cells to modify
            nbCellType = 0
            do iCell = 1, nbCell
                cellNumeIn = listCellNume(iCell)
                cellTypeIn = meshTypmailIn(cellNumeIn)
                if (cellTypeIn .eq. cellTypeToModify) then
                    nbCellType  = nbCellType + 1
                endif
            end do
            ASSERT(nbCellType .le. nbCell)

            do iCell = 1, nbCell
! ------------- Current cell
                cellNumeIn = listCellNume(iCell)
                cellTypeIn = meshTypmailIn(cellNumeIn)

! ------------- This type has to been modified => one node added, one cell modify
                if (cellTypeIn .eq. cellTypeToModify) then
                    nbCellModi  = nbCellModi + 1

                    do iCellModi = iad, iad + nbCellType - 1
                        ASSERT(iCellModi .le. nbCellIn)
                        creaNodeName(iCellModi) = prefNodeName
                    end do

                    ASSERT(iad .le. nbCellIn)
                    creaNodeNume(iad) = prefNodeNume

! ----------------- Save
                    ASSERT(nbCellModi .le. nbCellIn)
                    modiCellNume(nbCellModi) = cellNumeIn
                    modiCellType(nbCellModi) = cellTypeIn
                endif
            end do

            iad = iad + nbCellType
            call jedetr(jvCellNume)

! --------- No cells of this type to modify
            if (nbCellType .le. 0) then
                call utmess('A', 'MESH2_4', sk=option)
            endif
        end do
!
        if (nbCellModi .eq. 0) then
            call utmess('A', 'MESH2_5')
        endif

! ----- Each modified cell add ONE node
        nbNodeCrea = nbCellModi
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "CREA_MAILLE"
!
! --------------------------------------------------------------------------------------------------
!
    nbCellCrea = 0
    if (nbOccCreaMaille .ne. 0) then
        keywfact = 'CREA_MAILLE'
        call wkvect(crgrnu, 'V V I', nbCellIn, vi = creaCellNume)
        call wkvect(crgrno, 'V V K8', nbCellIn, vk8 = creaCellName)
        listCreaLength = nbCellIn
        nbCellCrea     = 0
        AS_ALLOCATE(vi = creaCellOccNb, size = nbOccCreaMaille)
        AS_ALLOCATE(vk24 = creaCellOccGrName, size = nbOccCreaMaille)

        do iocc = 1, nbOccCreaMaille
! --------- Get options from user for name of new cells
            call getCellOptionForName(keywfact     , iocc,&
                                      lPrefCellName, lPrefCellNume,&
                                      prefCellName , prefCellNume)

! --------- Get list of cells
            call getelem(meshIn, keywfact, iocc, 'F', jvCellNume, nbCell)
            call jeveuo(jvCellNume, 'L', vi = listCellNume)

! --------- Name of group of cells
            call getvtx(keywfact, 'NOM', iocc=iocc, scal=grCellName, nbret = n1)
            ASSERT(n1 .eq. 1)
            creaCellOccNb(iocc)     = nbCell
            creaCellOccGrName(iocc) = grCellName

            do iCell = 1, nbCell
! ------------- Cell to copy
                cellNumeIn = listCellNume(iCell)
                call jenuno(jexnum(meshIn//'.NOMMAI', cellNumeIn), cellNameIn)

! ------------- Create name for new cell
                cellNameOut = cellNameIn
                call createNameOfCell(cellNameOut  ,&
                                      lPrefCellName, lPrefCellNume,&
                                      prefCellName , prefCellNume)

! ------------- Check if name of cell exist
                call jenonu(jexnom(meshIn//'.NOMMAI', cellNameOut), cellNumeOut)
                if (cellNumeOut .ne. 0) then
                    call utmess('F', 'MESH2_2', sk = cellNameOut)
                endif

! ------------- A new cell
                nbCellCrea = nbCellCrea + 1

! ------------- Extend size of objects 
                if (nbCellCrea .gt. listCreaLength) then
                    call juveca(crgrno, 2*nbCellCrea)
                    call juveca(crgrnu, 2*nbCellCrea)
                    call jeveuo(crgrnu, 'E', vi = creaCellNume)
                    call jeveuo(crgrno, 'E', vk8 = creaCellName)
                    call jelira(crgrno, 'LONMAX', listCreaLength)
                endif

! ------------- Save name and index of new cell
                creaCellNume(nbCellCrea) = cellNumeIn
                creaCellName(nbCellCrea) = cellNameOut

            end do
            call jedetr(jvCellNume)
        end do
    endif
!
! - Each CREA_MAILLE create a group of cells
!
    nbGrCellFromCreaCell = nbOccCreaMaille
!
! --------------------------------------------------------------------------------------------------
!
!   For "CREA_POI1"
!
! --------------------------------------------------------------------------------------------------
!
    nbCellAddPoi1        = 0
    nbGrCellFromCreaPoi1 = 0
    cellShift            = 0
    if (nbOccCreaPoi1 .ne. 0) then
        keywfact = 'CREA_POI1'
        AS_ALLOCATE(vi = creaGrPoi1NbCell, size = nbOccCreaPoi1)
        AS_ALLOCATE(vk24 = creaGrPoi1GrName, size = nbOccCreaPoi1)
        AS_ALLOCATE(vk8 = creaGrPoi1CellName, size = nbOccCreaPoi1*nbNodeIn)
        AS_ALLOCATE(vl = creaPoi1Flag, size = nbNodeIn)
        AS_ALLOCATE(vk8 = creaPoi1Name, size = nbNodeIn)
        do iocc = 1, nbOccCreaPoi1

! --------- Get list of nodes to create POI1
            call getnode(meshIn, keywfact, iocc, 'F', jvNodeNume, nbNode)
            call jeveuo(jvNodeNume, 'L', vi = listNodeNume)
            do iNode = 1, nbNode
                nodeNumeIn = listNodeNume(iNode)
                if (.not. creaPoi1Flag(nodeNumeIn)) then
                    creaPoi1Flag(nodeNumeIn) = ASTER_TRUE
                endif
            end do

! --------- Get list of nodes to create GROUP_MA for POI1
            call getvtx(keywfact, 'NOM_GROUP_MA', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvtx(keywfact, 'NOM_GROUP_MA', iocc=iocc, scal = grCellName)
                creaGrPoi1GrName(iOcc) = grCellName
                nbGrCellFromCreaPoi1   = nbGrCellFromCreaPoi1 + 1
                creaGrPoi1NbCell(iOcc) = nbNode
                do iNode = 1, nbNode
                    nodeNumeIn = listNodeNume(iNode)
                    call jenuno(jexnum(meshIn//'.NOMNOE', nodeNumeIn), nodeNameIn)
                    cellNameOut = nodeNameIn
                    creaGrPoi1CellName(cellShift+iNode) = cellNameOut
                end do
                cellShift = cellShift + nbNode
            endif
            call jedetr(jvNodeNume)
        end do
! ----- Prepare objects for POI1 to add
        do iNode = 1, nbNodeIn
            nodeNumeIn = iNode
            if (.not. creaPoi1Flag(nodeNumeIn)) then
                cycle
            else
                nbCellAddPoi1 = nbCellAddPoi1 + 1
            endif
            call jenuno(jexnum(meshIn//'.NOMNOE', iNode), nodeNameIn)
            cellNameIn = nodeNameIn
            call jenonu(jexnom(meshIn//'.NOMMAI', cellNameIn), cellNumeIn)
            if (cellNumeIn .eq. 0) then
                creaPoi1Name(nbCellAddPoi1) = cellNameIn
            else
                call utmess('F', 'MESH1_14', sk = cellNameIn)
            endif
        end do
        AS_DEALLOCATE(vl = creaPoi1Flag)
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   Prepare base objects of mesh
!
! --------------------------------------------------------------------------------------------------
!
    call jedupo(meshIn//'.DIME'           , 'G', meshOut//'.DIME'           , ASTER_FALSE)
    call jedupo(meshIn//'.COORDO    .DESC', 'G', meshOut//'.COORDO    .DESC', ASTER_FALSE)
    call jedupo(meshIn//'.COORDO    .REFE', 'G', meshOut//'.COORDO    .REFE', ASTER_FALSE)
    call jedupo(meshIn//'.NOMACR'         , 'G', meshOut//'.NOMACR'         , ASTER_FALSE)
    call jedupo(meshIn//'.PARA_R'         , 'G', meshOut//'.PARA_R'         , ASTER_FALSE)
    call jedupo(meshIn//'.SUPMAIL'        , 'G', meshOut//'.SUPMAIL'        , ASTER_FALSE)
    call jedupo(meshIn//'.TYPL'           , 'G', meshOut//'.TYPL'           , ASTER_FALSE)
    call jedupo(meshIn//'.ABSC_CURV'      , 'G', meshOut//'.ABSC_CURV'      , ASTER_FALSE)
    call jedupo(meshIn//'.NOMACR', 'G', meshOut//'.NOMACR', ASTER_FALSE)
    call jedupo(meshIn//'.PARA_R', 'G', meshOut//'.PARA_R', ASTER_FALSE)
    call jedupo(meshIn//'.SUPMAIL', 'G', meshOut//'.SUPMAIL', ASTER_FALSE)
    call jedupo(meshIn//'.TYPL', 'G', meshOut//'.TYPL', ASTER_FALSE)
    call jedupo(meshIn//'.ABSC_CURV', 'G', meshOut//'.ABSC_CURV', ASTER_FALSE)

! - New sizes

    nbNodeOut   = nbNodeIn + nbNodeCrea
    nbCellOut   = nbCellIn + nbCellCrea + nbCellAddPoi1
    nbGrCellOut = nbGrCellIn + nbGrCellFromCreaCell + nbGrCellFromCreaPoi1
    nbGrNodeOut = nbGrNodeIn
    call jedupo(meshIn//'.DIME', 'G', meshOut//'.DIME', ASTER_FALSE)
    call jeveuo(meshOut//'.DIME', 'E', vi = meshDimeOut)
    meshDimeOut(1) = nbNodeOut
    meshDimeOut(3) = nbCellOut

! - Repertory of name of nodes
    if (nbNodeOut .eq. nbNodeIn) then
        call jedupo(meshIn//'.NOMNOE', 'G', meshOut//'.NOMNOE', ASTER_FALSE)
    else
        call jecreo(meshOut//'.NOMNOE', 'G N K8')
        call jeecra(meshOut//'.NOMNOE', 'NOMMAX', nbNodeOut, ' ')
! ----- Copy previous nodes
        do iNode = 1, nbNodeIn
            call jenuno(jexnum(meshIn//'.NOMNOE', iNode), nodeNameIn)
            nodeNameOut = nodeNameIn
            call jeexin(jexnom(meshOut//'.NOMNOE', nodeNameOut), iret)
            ASSERT(iret .eq. 0)
            call jecroc(jexnom(meshOut//'.NOMNOE', nodeNameOut))
        end do
! ----- Add new nodes
        do iNode = nbNodeIn+1, nbNodeOut
            call codent(creaNodeNume(iNode-nbNodeIn), 'G', knume)
            if (creaNodeName(iNode-nbNodeIn) .eq. creaNodeName(iNode-nbNodeIn+1)) then
                creaNodeNume(iNode-nbNodeIn+1) = creaNodeNume(iNode-nbNodeIn) + 1
            endif
            prefNodeName = creaNodeName(iNode-nbNodeIn)
            if (lxlgut(knume) + lxlgut(prefNodeName) .gt. 8) then
                call utmess('F', 'MESH2_1')
            endif
            nodeNameOut = prefNodeName(1:lxlgut(prefNodeName))//knume
            call jeexin(jexnom(meshOut//'.NOMNOE', nodeNameOut), nodeNumeOut)
            if (nodeNumeOut .eq. 0) then
                call jecroc(jexnom(meshOut//'.NOMNOE', nodeNameOut))
            else
                call utmess('F', 'MESH2_3', sk = nodeNameOut)
            endif
        end do
    endif

! - Coordinates of nodes
    call jedupo(meshIn//'.COORDO    .DESC', 'G', meshOut//'.COORDO    .DESC', ASTER_FALSE)
    call jedupo(meshIn//'.COORDO    .REFE', 'G', meshOut//'.COORDO    .REFE', ASTER_FALSE)
    call jeveuo(meshOut//'.COORDO    .REFE', 'E', vk24 = meshRefeOut)
    meshRefeOut(1) = meshOut
    if (nbNodeOut .eq. nbNodeIn) then
        call jedupo(meshIn//'.COORDO    .VALE', 'G', meshOut//'.COORDO    .VALE', ASTER_FALSE)
    else
        call jeveuo(meshIn//'.COORDO    .VALE', 'L', vr = meshValeIn)
        call wkvect(meshOut//'.COORDO    .VALE', 'G V R8', 3*nbNodeOut, vr = meshValeOut)
        meshValeOut(1:3*nbNodeIn) = meshValeIn(1:3*nbNodeIn)
        call jelira(meshIn//'.COORDO    .VALE', 'DOCU', cval=cdim)
        call jeecra(meshOut//'.COORDO    .VALE', 'DOCU', cval=cdim)
    endif
!
! ==================================================================================================
!
!   Modification of cells
!
! ==================================================================================================
!

! - Create cells
    call jecreo(meshOut//'.NOMMAI', 'G N K8')
    call jeecra(meshOut//'.NOMMAI', 'NOMMAX', nbCellOut, ' ')
    call wkvect(meshOut//'.TYPMAIL', 'G V I', nbCellOut, vi = meshTypmailOut)
    call jecrec(meshOut//'.CONNEX', 'G V I', 'NU', 'CONTIG', 'VARIABLE', nbCellOut)

! - Object to save parameters of cells
    AS_ALLOCATE(vi=nbNodeByCellOut, size=nbCellOut)
    AS_ALLOCATE(vi=connexAdr, size=nbCellOut)
    AS_ALLOCATE(vi=allCellNume, size=nbCellOut)

! - Copy previous cells
    connexLength = 0
    cellShift    = 0
    do iCell = 1, nbCellIn

! ----- Copy name of cell
        call jenuno(jexnum(meshIn//'.NOMMAI', iCell), cellNameIn)
        call jenonu(jexnom(meshIn//'.NOMMAI', cellNameIn), cellNumeIn)
        cellNameOut = cellNameIn
        call jecroc(jexnom(meshOut//'.NOMMAI', cellNameOut))
        call jenonu(jexnom(meshOut//'.NOMMAI', cellNameOut), cellNumeOut)

! ----- Copy type of cell
        meshTypmailOut(cellNumeOut) = meshTypmailIn(cellNumeIn)

! ----- Copy connexity of cell
        call jelira(jexnum(meshIn//'.CONNEX', cellNumeIn), 'LONMAX', nbNodeInCellIn)
        call jeveuo(jexnum(meshIn//'.CONNEX', cellNumeIn), 'L', jvConnexIn)
        nbNodeInCellOut = nbNodeInCellIn
        do iNode = 1, nbNodeCrea
            if (iCell .eq. modiCellNume(iNode)) then
                nbNodeInCellOut = nbNodeInCellIn + 1
                goto 160
            endif
        end do
160     continue
        connexLength = connexLength + nbNodeInCellOut

! ----- Save new parameters of cell
        nbNodeByCellOut(iCell) = nbNodeInCellOut
        connexAdr(iCell)       = jvConnexIn
        allCellNume(iCell)     = cellNumeOut
    end do
    cellShift = cellShift + nbCellIn

! - Create cell from "CREA_MAILLE"
    do iCell = 1, nbCellCrea
        cellNameOut = creaCellName(iCell)
        cellNumeIn  = creaCellNume(iCell)

! ----- Create name of new cell
        call jeexin(jexnom(meshOut//'.NOMMAI', cellNameOut), cellNumeOut)
        if (cellNumeOut .eq. 0) then
            call jecroc(jexnom(meshOut//'.NOMMAI', cellNameOut))
        else
            call utmess('F', 'MESH2_2', sk = cellNameOut)
        endif
        call jenonu(jexnom(meshOut//'.NOMMAI', cellNameOut), cellNumeOut)

! ----- Copy type of new cell
        meshTypmailOut(cellNumeOut) = meshTypmailIn(cellNumeIn)

! ----- Get connexity of old cell
        call jelira(jexnum(meshIn//'.CONNEX', cellNumeIn), 'LONMAX', nbNodeInCellIn)
        call jeveuo(jexnum(meshIn//'.CONNEX', cellNumeIn), 'L', jvConnexIn)
        connexLength    = connexLength + nbNodeInCellIn
        nbNodeInCellOut = nbNodeInCellIn

! ----- Save new parameters of cell
        nbNodeByCellOut(cellShift+iCell) = nbNodeInCellOut
        connexAdr(cellShift+iCell)       = jvConnexIn
        allCellNume(cellShift+iCell)     = cellNumeOut

    end do

! - Set new connexity of cells
    connexLength = connexLength + nbCellAddPoi1
    call jeecra(meshOut//'.CONNEX', 'LONT', connexLength)

! - Copy connexity from previous cells
    cellShift = 0
    do iCell = 1, nbCellIn
        nbNodeInCellOut = nbNodeByCellOut(cellShift+iCell)
        jvConnexIn      = connexAdr(cellShift+iCell)
        cellNumeOut     = allCellNume(cellShift+iCell)
        call jeecra(jexnum(meshOut//'.CONNEX', cellNumeOut), 'LONMAX', nbNodeInCellOut)
        call jeveuo(jexnum(meshOut//'.CONNEX', cellNumeOut), 'E', jvConnexOut)
        do iNode = 1, nbNodeInCellOut
            zi(jvConnexOut - 1 + iNode) = zi(jvConnexIn - 1 + iNode)
        end do
    end do
    cellShift = cellShift + nbCellIn

! - Create connexity for CREA_MAILLE
    do iCell = 1, nbCellCrea
        nbNodeInCellOut = nbNodeByCellOut(cellShift+iCell)
        jvConnexIn      = connexAdr(cellShift+iCell)
        cellNumeOut     = allCellNume(cellShift+iCell)
        call jeecra(jexnum(meshOut//'.CONNEX', cellNumeOut), 'LONMAX', nbNodeInCellOut)
        call jeveuo(jexnum(meshOut//'.CONNEX', cellNumeOut), 'E', jvConnexOut)
        do iNode = 1, nbNodeInCellOut
            zi(jvConnexOut - 1 + iNode) = zi(jvConnexIn - 1 + iNode)
        end do
    end do

! - Create cell for CREA_POI1
    do iCell = 1, nbCellAddPoi1
        cellNameOut = creaPoi1Name(iCell)
        nodeNameOut = cellNameOut

! ----- Create name of new cell
        call jeexin(jexnom(meshOut//'.NOMMAI', cellNameOut), cellNumeOut)
        if (cellNumeOut .eq. 0) then
            call jecroc(jexnom(meshOut//'.NOMMAI', cellNameOut))
        else
            call utmess('F', 'MESH2_2', sk = cellNameOut)
        endif
        call jenonu(jexnom(meshOut//'.NOMMAI', cellNameOut), cellNumeOut)

! ----- Set type of cell
        meshTypmailOut(cellNumeOut) = MT_POI1

! ----- Set connectivity
        nbNodeInCellOut = 1
        call jeecra(jexnum(meshOut//'.CONNEX', cellNumeOut), 'LONMAX', nbNodeInCellOut)
        call jeveuo(jexnum(meshOut//'.CONNEX', cellNumeOut), 'E', jvConnexOut)
        call jenonu(jexnom(meshOut//'.NOMNOE', nodeNameOut), nodeNumeOut)
        zi(jvConnexOut) = nodeNumeOut

    end do
    AS_DEALLOCATE(vi=nbNodeByCellOut)
    AS_DEALLOCATE(vi=connexAdr)
    AS_DEALLOCATE(vi=allCellNume)
!
! ==================================================================================================
!
!   Modification of groups of cells
!
! ==================================================================================================
!

! - Create repertory of names for groups of cells
    if (nbGrCellOut .ne. 0) then
        call jecreo(meshOut//'.PTRNOMMAI', 'G N K24')
        call jeecra(meshOut//'.PTRNOMMAI', 'NOMMAX', nbGrCellOut, ' ')
        call jecrec(meshOut//'.GROUPEMA', 'G V I', 'NO '//meshOut//'.PTRNOMMAI',&
                   'DISPERSE', 'VARIABLE', nbGrCellOut)
    endif

! - Copy previous groups of cells
    if (nbGrCellOut .ne. 0) then
        do iGrCell = 1, nbGrCellIn

! --------- Get group from input mesh
            call jenuno(jexnum(meshIn//'.GROUPEMA', iGrCell), grCellName)
            call jeveuo(jexnum(meshIn//'.GROUPEMA', iGrCell), 'L', vi = cellInGrIn)
            call jelira(jexnum(meshIn//'.GROUPEMA', iGrCell), 'LONMAX', nbCellInGrIn)
            call jelira(jexnum(meshIn//'.GROUPEMA', iGrCell), 'LONUTI', nbCellInGrIn)

! --------- Create group in output mesh
            call jeexin(jexnom(meshOut//'.GROUPEMA', grCellName), iret)
            ASSERT(iret .eq. 0)
            call jecroc(jexnom(meshOut//'.GROUPEMA', grCellName))
            nbCellInGrOut = nbCellInGrIn
            call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONMAX', max(nbCellInGrOut, 1))
            call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONUTI', nbCellInGrOut)
            call jeveuo(jexnom(meshOut//'.GROUPEMA', grCellName), 'E', vi = cellInGrOut)
            cellInGrOut(1:nbCellInGrOut) = cellInGrIn(1:nbCellInGrIn)
        end do
    endif

! - Create new groups of cells for CREA_MAILLE
    cellShift = 0
    do iOcc = 1, nbOccCreaMaille
        grCellName    = creaCellOccGrName(iOcc)
        nbCellInGrOut = creaCellOccNb(iOcc)

! ----- Create name of group in repertory
        call jeexin(jexnom(meshOut//'.GROUPEMA', grCellName), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(meshOut//'.GROUPEMA', grCellName))
        else
            call utmess('F', 'MESH1_20', sk = grCellName)
        endif

! ----- Create group in output mesh
        call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONMAX', max(nbCellInGrOut, 1))
        call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONUTI', nbCellInGrOut)
        call jeveuo(jexnom(meshOut//'.GROUPEMA', grCellName), 'E', vi = cellInGrOut)
        do iCell = 1, nbCellInGrOut
            cellNameOut = creaCellName(iCell+cellShift)
            call jenonu(jexnom(meshOut//'.NOMMAI', cellNameOut), cellInGrOut(iCell))
        end do
        cellShift = cellShift + nbCellInGrOut
    end do
    ASSERT(cellShift .eq. nbCellCrea)

! - Create new groups of cells for CREA_POI1
    cellShift = 0
    do iOcc = 1, nbOccCreaPoi1
        grCellName    = creaGrPoi1GrName(iOcc)
        nbCellInGrOut = creaGrPoi1NbCell(iOcc)
        if (grCellName .ne. ' ') then

! --------- Create name of group in repertory
            call jeexin(jexnom(meshOut//'.GROUPEMA', grCellName), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(meshOut//'.GROUPEMA', grCellName))
            else
                call utmess('F', 'MESH1_20', sk = grCellName)
            endif

! --------- Create group in output mesh
            call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONMAX', max(nbCellInGrOut, 1))
            call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONUTI', nbCellInGrOut)
            call jeveuo(jexnom(meshOut//'.GROUPEMA', grCellName), 'E', vi = cellInGrOut)
            do iCell = 1, nbCellInGrOut
                cellNameOut = creaGrPoi1CellName(iCell+cellShift)
                call jenonu(jexnom(meshOut//'.NOMMAI', cellNameOut), cellInGrOut(iCell))
            end do
            cellShift = cellShift + nbCellInGrOut
        endif
    end do
!
! ==================================================================================================
!
!   Modification of groups of nodes
!
! ==================================================================================================
!
    if (nbGrNodeOut .ne. 0) then
        call jecreo(meshOut//'.PTRNOMNOE', 'G N K24')
        call jeecra(meshOut//'.PTRNOMNOE', 'NOMMAX', nbGrNodeOut, ' ')
        call jecrec(meshOut//'.GROUPENO', 'G V I', 'NO '//meshOut//'.PTRNOMNOE',&
                    'DISPERSE', 'VARIABLE', nbGrNodeOut)
        do iGrNode = 1, nbGrNodeIn
! --------- Get group from input mesh
            call jenuno(jexnum(meshIn//'.GROUPENO', iGrNode), grNodeName)
            call jeveuo(jexnum(meshIn//'.GROUPENO', iGrNode), 'L', vi = nodeInGrIn)
            call jelira(jexnum(meshIn//'.GROUPENO', iGrNode), 'LONUTI', nbNodeInGrIn)

! --------- Create group in output mesh
            nbNodeInGrOut = nbNodeInGrIn
            call jeexin(jexnom(meshOut//'.GROUPENO', grNodeName), iret)
            ASSERT(iret .eq. 0)
            call jecroc(jexnom(meshOut//'.GROUPENO', grNodeName))
            call jeecra(jexnom(meshOut//'.GROUPENO', grNodeName), 'LONMAX', max(nbNodeInGrOut, 1))
            call jeecra(jexnom(meshOut//'.GROUPENO', grNodeName), 'LONUTI', nbNodeInGrOut)
            call jeveuo(jexnom(meshOut//'.GROUPENO', grNodeName), 'E', vi = nodeInGrOut)
            nodeInGrOut(1:nbNodeInGrOut) = nodeInGrIn(1:nbNodeInGrIn)
        end do
    endif
!
! - Add new nodes from modified cells
!
    if (nbCellModi .ne. 0) then
        call cmmoma(meshOut, nbCellModi, modiCellNume, modiCellType, nbNodeIn)
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "DECOUPE_LAC"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccDecoupeLac .gt. 0) then
        ASSERT(nbOccDecoupeLac.eq.1)
        call cm_dclac(meshIn, meshOut)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "REPERE"
!
! --------------------------------------------------------------------------------------------------
!
    if (nbOccRepere .ne. 0) then
        ASSERT(nbOccRepere .eq. 1)
        keywfact = 'REPERE'
        call getvid(keywfact, 'TABLE', iocc=1, nbval=0, nbret=ntab)
        ntab = abs(ntab)
        ASSERT(ntab .eq. 1)
        call getvid(keywfact, 'TABLE', iocc=1, scal=table)
        call getvtx(keywfact, 'NOM_ORIG', iocc=1, nbval=0, nbret=nori)
        if (nori .ne. 0) then
            call getvtx(keywfact, 'NOM_ORIG', iocc=1, scal=nomori)
            if (nomori .eq. 'CDG') then
                call chcoma(table, meshOut)
            else if (nomori.eq.'TORSION') then
                call chcomb(table, meshOut)
            else
                ASSERT(ASTER_FALSE)
            endif
        endif
        goto 350
    endif
!
350 continue

! - Add title in mesh datastructure
    call titre()

! - Update parameters for modified mesh (bounding box and dimensions)
    call cargeo(meshOut)

! - Verbose
    call infoma(meshOut)

! - Clean
    AS_DEALLOCATE(vi = modiCellNume)
    AS_DEALLOCATE(vi = modiCellType)
    AS_DEALLOCATE(vk8 = creaNodeName)
    AS_DEALLOCATE(vi = creaNodeNume)
    AS_DEALLOCATE(vi = creaCellOccNb)
    AS_DEALLOCATE(vk24 = creaCellOccGrName)
    AS_DEALLOCATE(vi = creaGrPoi1NbCell)
    AS_DEALLOCATE(vk24 = creaGrPoi1GrName)
    AS_DEALLOCATE(vk8 = creaGrPoi1CellName)
    AS_DEALLOCATE(vk8 = creaPoi1Name)
!
    call jedema()
!
end subroutine
