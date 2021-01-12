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
#include "asterfort/cpclma.h"
#include "asterfort/dismoi.h"
#include "asterfort/eclpgm.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
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
    integer :: i, nbmc, iret, nbma, iqtr
    integer :: n1
    integer :: n1b, jlima, nbmoma
    parameter(nbmc=5)
    real(kind=8) :: epais
    character(len=4) :: answer
    character(len=4) :: cdim
    character(len=8) :: meshIn, meshOut, model, geofi
    character(len=8) :: nomori, plan, trans, knume
    character(len=16) :: option, keywfact
    character(len=16) :: kbi1, kbi2
    character(len=16) :: motfac, tymocl(nbmc), motcle(nbmc)
    character(len=19) :: table, ligrel
    character(len=19), parameter :: k19void = ' '
    integer :: nbMeshIn
    character(len=24) :: nommai, grpmai, typmai, connex, nodime, grpnoe, nomnoe
    character(len=24) :: cooval, coodsc, cooref, nomjv
    character(len=24) :: nommav, grpmav, typmav, connev, nodimv, grpnov, nomnov
    character(len=24) :: coovav, coodsv, coorev
    character(len=24) :: crgrnu, crgrno
    character(len=24) :: grCellName, grNodeName, valk(2), nogma, gpptnm, gpptnn
    integer :: iagma, iatyma, ino, j
    integer :: jgg, jmail
    integer :: jnoeu, jnono, jnpt, jopt, jtom, jtrno, jvale, jvg, kvale
    integer :: nbgma, nbgrmt
    integer :: nbgrno, nbmain, nbmaj3, nbno
    integer :: nbpt, nbptt, nori, ntab, ntpoi
    integer :: ibid, ifm, jdime
    integer :: jrefe, jtypmv
    integer :: nbmaiv, niv, k, jgeofi
    integer :: dimcon, decala
    integer :: cellNume, creaCellNume, cellType, cellTypeToModify
    integer :: nodeNume
    integer :: iCell, iNode, iCellModi, iGrCell
    integer :: shiftCell
    integer :: iocc, nbOcc
    character(len=24), parameter :: jvCellNume = '&&OP0167.LISTCELL'
    integer :: nbCell, nbCellCrea, nbCellModi, nbCellType
    integer :: nbGrCellIn, nbGrCellOut
    integer :: nbNodeCrea, nbNodeIn, nbNodeOut
    integer :: nbField
    integer :: nbOccDecoupeLac, nbOccEclaPg, nbGeomFibre, nbOccCreaFiss, nbOccLineQuad
    integer :: nbOccQuadLine, nbOccModiMaille, nbOccCoquVolu, nbOccRestreint, nbOccRepere
    integer :: iOccQuadTria, iad
    integer :: nbOccCreaPoi1, nbOccCreaMaille
    real(kind=8) :: shrink, lonmin
    aster_logical :: lpb
    character(len=8) :: cellName, nodeName, creaCellName
    aster_logical :: lPrefCellName, lPrefCellNume, lPrefNodeName, lPrefNodeNume
    integer :: prefCellNume, prefNodeNume, prefNume
    character(len=8) :: prefCellName, prefNodeName
    integer, pointer :: modiCellNume(:) => null()
    integer, pointer :: modiCellType(:) => null()
    character(len=8), pointer :: addNodeName(:) => null()
    integer, pointer :: addNodeNume(:) => null()
    integer, pointer :: listCellNume(:) => null()
    character(len=16), pointer :: listField(:) => null()
    integer, pointer :: adrjvx(:) => null()
    integer, pointer :: nbnoma(:) => null()
    integer, pointer :: nbnomb(:) => null()
    integer, pointer :: nomnum(:) => null()
    integer, pointer :: listCreaNume(:) => null()
    character(len=8), pointer :: listCreaName(:) => null()
    integer, pointer :: listCreaOccNb(:) => null()
    character(len=24), pointer :: listCreaOccGrName(:) => null()
    integer :: listCreaLength
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
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
!
! - Main datastructure
!
    call getres(meshOut, kbi1, kbi2)
    call getvid(' ', 'MAILLAGE', scal = meshIn, nbret = nbMeshIn)
    if (nbMeshIn .ne. 0) then
        if (isParallelMesh(meshIn)) then
            call utmess('F', 'MESH1_22')
        end if
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
        call jeveuo(geofi//'.GFMA', 'L', jgeofi)
        call copisd('MAILLAGE', 'G', zk8(jgeofi), meshOut)
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
        if (nbCell .ne. nbmaiv) then
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
            if (nbCell .ne. nbmaiv) then
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
        if (nbCell .ne. nbmaiv) then
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
    call getfac('MODI_HHO', nbmoma)
    if (nbmoma .gt. 0) then
        ASSERT(nbmoma.eq.1)
!
        call getvtx('MODI_HHO', 'GROUP_MA', iocc=1, nbval=0, nbret=n1b)
!
        call getvtx("MODI_HHO", 'PREF_NOEUD', iocc=1, scal=prefNodeName)
        call getvis("MODI_HHO", 'PREF_NUME', iocc=1, scal=prefNodeNume)
!
        motcle(1)='GROUP_MA'
        motcle(2)='TOUT'
        nomjv='&&OP0167.LISTE_MA'
        call reliem(' ', meshIn, 'NU_MAILLE', 'MODI_HHO', 1,&
                    2, motcle, motcle, nomjv, nbma)

        if (nbma .ne. nbmaiv) then
            call utmess('A', 'MESH1_4', sk='MODI_HHO')
        endif

        call jeveuo(nomjv, 'L', jlima)
        call jeexin(meshIn//'.NOMACR', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_7')
        endif
        call jeexin(meshIn//'.ABSC_CURV', iret)
        if (iret .ne. 0) then
            call utmess('F', 'MESH1_8')
        endif
!
        call cmhho(meshIn, meshOut, nbma, zi(jlima), prefNodeName, prefNodeNume)
!
        goto 350
!
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
            if (nbCell .ne. nbmaiv) then
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
!   Copy base objects of mesh
!
! --------------------------------------------------------------------------------------------------
!
    nommav=meshIn//'.NOMMAI         '
    nomnov=meshIn//'.NOMNOE         '
    typmav=meshIn//'.TYPMAIL        '
    connev=meshIn//'.CONNEX         '
    grpmav=meshIn//'.GROUPEMA       '
    grpnov=meshIn//'.GROUPENO       '
    nodimv=meshIn//'.DIME           '
    coovav=meshIn//'.COORDO    .VALE'
    coodsv=meshIn//'.COORDO    .DESC'
    coorev=meshIn//'.COORDO    .REFE'
!
    nommai=meshOut//'.NOMMAI         '
    nomnoe=meshOut//'.NOMNOE         '
    typmai=meshOut//'.TYPMAIL        '
    connex=meshOut//'.CONNEX         '
    grpmai=meshOut//'.GROUPEMA       '
    grpnoe=meshOut//'.GROUPENO       '
    nodime=meshOut//'.DIME           '
    cooval=meshOut//'.COORDO    .VALE'
    coodsc=meshOut//'.COORDO    .DESC'
    cooref=meshOut//'.COORDO    .REFE'
    gpptnm=meshOut//'.PTRNOMMAI'
    gpptnn=meshOut//'.PTRNOMNOE'
!
!
    call jedupo(nodimv, 'G', nodime, .false._1)
    call jedupo(coodsv, 'G', coodsc, .false._1)
    call jedupo(coorev, 'G', cooref, .false._1)
    call jedupo(meshIn//'.NOMACR', 'G', meshOut//'.NOMACR', .false._1)
    call jedupo(meshIn//'.PARA_R', 'G', meshOut//'.PARA_R', .false._1)
    call jedupo(meshIn//'.SUPMAIL', 'G', meshOut//'.SUPMAIL', .false._1)
    call jedupo(meshIn//'.TYPL', 'G', meshOut//'.TYPL', .false._1)
    call jedupo(meshIn//'.ABSC_CURV', 'G', meshOut//'.ABSC_CURV', .false._1)
!
    call jeveuo(cooref, 'E', jrefe)
    zk24(jrefe)=meshOut
!
    call jeveuo(nodime, 'E', jdime)
    nbNodeIn = zi(jdime)
    nbmaiv   = zi(jdime+3-1)
!
    call jeveuo(typmav, 'L', jtypmv)
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
        AS_ALLOCATE(vi = modiCellNume, size = nbmaiv)
        AS_ALLOCATE(vi = modiCellType, size = nbmaiv)
        AS_ALLOCATE(vk8 = addNodeName, size = nbmaiv)
        AS_ALLOCATE(vi  = addNodeNume, size = nbmaiv)

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
                cellNume = listCellNume(iCell)
                cellType = zi(jtypmv-1+cellNume)
                if (cellType .eq. cellTypeToModify) then
                    nbCellType  = nbCellType + 1
                endif
            end do
            ASSERT(nbCellType .le. nbCell)

            do iCell = 1, nbCell
! ------------- Current cell
                cellNume = listCellNume(iCell)
                cellType = zi(jtypmv-1+cellNume)

! ------------- This type has to been modified => one node added, one cell modify
                if (cellType .eq. cellTypeToModify) then
                    nbCellModi  = nbCellModi + 1

                    do iCellModi = iad, iad + nbCellType - 1
                        ASSERT(iCellModi .le. nbmaiv)
                        addNodeName(iCellModi) = prefNodeName
                    end do

                    ASSERT(iad .le. nbmaiv)
                    addNodeNume(iad) = prefNodeNume

! ----------------- Save
                    ASSERT(nbCellModi .le. nbmaiv)
                    modiCellNume(nbCellModi) = cellNume
                    modiCellType(nbCellModi) = cellType
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
        crgrnu   = '&&OP0167.CR_GR.NUM'
        crgrno   = '&&OP0167.CR_GR.NOM'
        call wkvect(crgrnu, 'V V I', nbmaiv, vi = listCreaNume)
        call wkvect(crgrno, 'V V K8', nbmaiv, vk8 = listCreaName)
        listCreaLength = nbmaiv
        nbCellCrea     = 0
        AS_ALLOCATE(vi = listCreaOccNb, size = nbOccCreaMaille)
        AS_ALLOCATE(vk24 = listCreaOccGrName, size = nbOccCreaMaille)

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
            listCreaOccNb(iocc)     = nbCell
            listCreaOccGrName(iocc) = grCellName

            do iCell = 1, nbCell
! ------------- Cell to copy
                cellNume = listCellNume(iCell)
                call jenuno(jexnum(meshIn//'.NOMMAI', cellNume), cellName)

! ------------- Create name for new cell
                creaCellName = cellName
                call createNameOfCell(creaCellName,&
                                      lPrefCellName, lPrefCellNume,&
                                      prefCellName , prefCellNume)

! ------------- Check if name of cell exist
                call jenonu(jexnom(meshIn//'.NOMMAI', creaCellName), creaCellNume)
                if (creaCellNume .ne. 0) then
                    call utmess('F', 'MESH2_2', sk = creaCellName)
                endif

! ------------- A new cell
                nbCellCrea = nbCellCrea + 1

! ------------- Extend size of objects 
                if (nbCellCrea .gt. listCreaLength) then
                    call juveca(crgrno, 2*nbCellCrea)
                    call juveca(crgrnu, 2*nbCellCrea)
                    call jeveuo(crgrnu, 'E', vi = listCreaNume)
                    call jeveuo(crgrno, 'E', vk8 = listCreaName)
                    call jelira(crgrno, 'LONMAX', listCreaLength)
                endif

! ------------- Save name and index of new cell
                listCreaNume(nbCellCrea) = cellNume
                listCreaName(nbCellCrea) = creaCellName

            end do
            call jedetr(jvCellNume)
        end do
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "CREA_POI1"
!
! --------------------------------------------------------------------------------------------------
!
    nbmaj3=0
    if (nbOccCreaPoi1 .ne. 0) then
        call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ntpoi)
!
!        -- RECUPERATION DE LA LISTE DES NOEUD :
        nomjv='&&OP0167.LISTE_NO'
        motfac='CREA_POI1'
        motcle(1)='NOEUD'
        tymocl(1)='NOEUD'
        motcle(2)='GROUP_NO'
        tymocl(2)='GROUP_NO'
        motcle(3)='MAILLE'
        tymocl(3)='MAILLE'
        motcle(4)='GROUP_MA'
        tymocl(4)='GROUP_MA'
        motcle(5)='TOUT'
        tymocl(5)='TOUT'
!
        call wkvect('&&OP0167.IND_NOEUD', 'V V I', nbNodeIn, jtrno)
        call wkvect('&&OP0167.NOM_NOEUD', 'V V K8', nbNodeIn, jnono)
!
        do iocc = 1, nbOccCreaPoi1
            call reliem(' ', meshIn, 'NO_NOEUD', motfac, iocc,&
                        nbmc, motcle, tymocl, nomjv, nbno)
            call jeveuo(nomjv, 'L', jnoeu)
            do i = 0, nbno-1
                call jenonu(jexnom(nomnov, zk8(jnoeu+i)), ino)
                zi(jtrno-1+ino)=1
            end do
        end do
!
!        --- VERIFICATION QUE LE NOM N'EXISTE PAS ET COMPTAGE---
        do iCell = 1, nbNodeIn
            if (zi(jtrno+iCell-1) .eq. 0) goto 110
            call jenuno(jexnum(nomnov, iCell), cellName)
            call jenonu(jexnom(nommav, cellName), ibid)
            if (ibid .eq. 0) then
                nbmaj3=nbmaj3+1
                zk8(jnono-1+nbmaj3)=cellName
            else
                valk(1)=cellName
                call utmess('F', 'ALGELINE4_43', nk=1, valk=valk)
            endif
110         continue
        end do
    endif
!
! ----------------------------------------------------------------------
!          ON AGRANDIT LE '.NOMNOE' ET LE '.COORDO    .VALE'
! ----------------------------------------------------------------------
!
    if (nbNodeCrea .ne. 0) then
        nbNodeOut=nbNodeIn+nbNodeCrea
        zi(jdime)=nbNodeOut
!
        call jecreo(nomnoe, 'G N K8')
        call jeecra(nomnoe, 'NOMMAX', nbNodeOut, ' ')
        do ino = 1, nbNodeIn
            call jenuno(jexnum(nomnov, ino), nodeName)
            call jeexin(jexnom(nomnoe, nodeName), iret)
            call jecroc(jexnom(nomnoe, nodeName))
        end do
        do ino = nbNodeIn+1, nbNodeOut
            call codent(addNodeNume(ino-nbNodeIn), 'G', knume)

            if (addNodeName(ino-nbNodeIn) .eq. addNodeName(ino-nbNodeIn+1)) then
                addNodeNume(ino-nbNodeIn+1) = addNodeNume(ino-nbNodeIn) + 1
            endif

            prefNodeName = addNodeName(ino-nbNodeIn)
            if (lxlgut(knume) + lxlgut(prefNodeName) .gt. 8) then
                call utmess('F', 'MESH2_1')
            endif

            nodeName = prefNodeName(1:lxlgut(prefNodeName))//knume
            call jeexin(jexnom(nomnoe, nodeName), nodeNume)
            if (nodeNume .eq. 0) then
                call jecroc(jexnom(nomnoe, nodeName))
            else
                call utmess('F', 'MESH2_3', sk = nodeName)

            endif
        end do

!
        call jeveuo(coovav, 'L', jvale)
        call wkvect(cooval, 'G V R8', 3*nbNodeOut, kvale)
        do i = 0, 3*nbNodeIn-1
            zr(kvale+i)=zr(jvale+i)
        end do
        call jelira(coovav, 'DOCU', cval=cdim)
        call jeecra(cooval, 'DOCU', cval=cdim)
    else
        call jedupo(nomnov, 'G', nomnoe, .false._1)
        call jedupo(coovav, 'G', cooval, .false._1)
    endif
!
! ----------------------------------------------------------------------
!         ON AGRANDIT LE '.NOMMAI' ET LE '.CONNEX'
! ----------------------------------------------------------------------
!
    nbmain=nbmaiv+nbCellCrea+nbmaj3
!
    zi(jdime+3-1)=nbmain
    call jecreo(nommai, 'G N K8')
    call jeecra(nommai, 'NOMMAX', nbmain, ' ')
!
    call wkvect(typmai, 'G V I', nbmain, iatyma)
!
    call jecrec(connex, 'G V I', 'NU', 'CONTIG', 'VARIABLE',&
                nbmain)
!
    AS_ALLOCATE(vi=nbnoma, size=nbmain)
    AS_ALLOCATE(vi=nbnomb, size=nbmaiv)
    AS_ALLOCATE(vi=adrjvx, size=nbmain)
    AS_ALLOCATE(vi=nomnum, size=nbmain)
    dimcon = 0
    decala = 0
    do iCell = 1, nbmaiv
        call jenuno(jexnum(nommav, iCell), cellName)
        call jeexin(jexnom(nommai, cellName), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(nommai, cellName))
        else
            call utmess('F', 'ALGELINE4_7', sk=cellName)
        endif
!
        call jenonu(jexnom(nommav, cellName), cellNume)
        jtom=jtypmv-1+cellNume
        call jenonu(jexnom(nommai, cellName), cellNume)
        zi(iatyma-1+cellNume)=zi(jtom)
!
        call jenonu(jexnom(nommav, cellName), cellNume)
        call jelira(jexnum(connev, cellNume), 'LONMAX', nbpt)
        call jeveuo(jexnum(connev, cellNume), 'L', jopt)
        nbptt=nbpt
        do iNode = 1, nbNodeCrea
            if (iCell .eq. modiCellNume(iNode)) then
                nbptt=nbpt+1
                goto 160
!
            endif
        end do
160     continue
        call jenonu(jexnom(nommai, cellName), cellNume)
        dimcon = dimcon+nbptt
        nbnoma(iCell) = nbptt
        nbnomb(iCell) = nbpt
        adrjvx(iCell) = jopt
        nomnum(iCell) = cellNume
    end do
!
    decala = decala + nbmaiv
!
    do iCell = 1, nbCellCrea
        creaCellName = listCreaName(iCell)
        cellNume     = listCreaNume(iCell)

! ----- Create name of new cell
        call jeexin(jexnom(nommai, creaCellName), creaCellNume)
        if (creaCellNume .eq. 0) then
            call jecroc(jexnom(nommai, creaCellName))
        else
            call utmess('F', 'ALGELINE4_7', sk=creaCellName)
        endif
        call jenonu(jexnom(nommai, creaCellName), creaCellNume)
        ASSERT(creaCellNume .gt. 0)

! ----- Copy type of new cell
        zi(iatyma-1+creaCellNume)=zi(jtypmv-1+cellNume)

! ----- Get connexity of old cell
        call jelira(jexnum(connev, cellNume), 'LONMAX', nbpt)
        call jeveuo(jexnum(connev, cellNume), 'L', jopt)
        dimcon = dimcon+nbpt
        nbnoma(1+decala+iCell-1) = nbpt
        adrjvx(1+decala+iCell-1) = jopt
        nomnum(1+decala+iCell-1) = creaCellNume
    end do
!
    dimcon = dimcon+nbmaj3
    call jeecra(connex, 'LONT', dimcon)
!
    decala = 0
    do iCell = 1, nbmaiv
        nbptt = nbnoma(1+decala+iCell-1)
        nbpt = nbnomb(1+decala+iCell-1)
        jopt = adrjvx(1+decala+iCell-1)
        ibid = nomnum(1+decala+iCell-1)
        call jeecra(jexnum(connex, ibid), 'LONMAX', nbptt)
        call jeveuo(jexnum(connex, ibid), 'E', jnpt)
        do ino = 0, nbpt-1
            zi(jnpt+ino)=zi(jopt+ino)
        end do
    end do
    decala = decala + nbmaiv

! - Add new cells for CREA_MAILLE
    do iCell = 1, nbCellCrea
        nbpt = nbnoma(1+decala+iCell-1)
        jopt = adrjvx(1+decala+iCell-1)
        creaCellNume = nomnum(1+decala+iCell-1)
        call jeecra(jexnum(connex, creaCellNume), 'LONMAX', nbpt)
        call jeveuo(jexnum(connex, creaCellNume), 'E', jnpt)
        do ino = 0, nbpt-1
            zi(jnpt+ino)=zi(jopt+ino)
        end do
    end do

!
    do iCell = 1, nbmaj3
        cellName=zk8(jnono+iCell-1)
        call jenonu(jexnom(nommai, cellName), ibid)
        if (ibid .ne. 0) goto 230
        call jeexin(jexnom(nommai, cellName), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(nommai, cellName))
        else
            valk(1)=cellName
            call utmess('F', 'ALGELINE4_7', sk=valk(1))
        endif
!
        call jenonu(jexnom(nommai, cellName), ibid)
        if (ibid .eq. 0) then
            call utmess('F', 'ALGELINE3_6', sk=cellName)
        endif
        zi(iatyma-1+ibid)=ntpoi
!
        call jeecra(jexnum(connex, ibid), 'LONMAX', 1)
        call jeveuo(jexnum(connex, ibid), 'E', jnpt)
        call jenonu(jexnom(nomnoe, cellName), zi(jnpt))
230     continue
    end do
    AS_DEALLOCATE(vi=nbnoma)
    AS_DEALLOCATE(vi=nbnomb)
    AS_DEALLOCATE(vi=adrjvx)
    AS_DEALLOCATE(vi=nomnum)
!
! ==================================================================================================
!
!   Modification of groups of cells
!
! ==================================================================================================
!
    call jeexin(meshIn//'.GROUPEMA', iret)
    if (iret .eq. 0) then
        nbGrCellIn = 0
    else
        call jelira(meshIn//'.GROUPEMA', 'NOMUTI', nbGrCellIn)
    endif
    nbGrCellOut = nbGrCellIn + nbOccCreaMaille

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
            call jenuno(jexnum(meshIn//'.GROUPEMA', iGrCell), grCellName)
            call jeexin(jexnom(meshOut//'.GROUPEMA', grCellName), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(meshOut//'.GROUPEMA', grCellName))
            else
                call utmess('F', 'ALGELINE4_9', sk=grCellName)
            endif
            call jeveuo(jexnum(meshIn//'.GROUPEMA', iGrCell), 'L', jvg)
            call jelira(jexnum(meshIn//'.GROUPEMA', iGrCell), 'LONMAX', nbCell)
            call jelira(jexnum(meshIn//'.GROUPEMA', iGrCell), 'LONUTI', nbCell)
            call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONMAX', max(nbCell, 1))
            call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONUTI', nbCell)
            call jeveuo(jexnom(meshOut//'.GROUPEMA', grCellName), 'E', jgg)
            do iCell = 1, nbCell
                zi(jgg-1+iCell)=zi(jvg-1+iCell)
            end do
        end do
    endif

! - Create new groups of cells for CREA_MAILLE
    shiftCell = 0
    do iOcc = 1, nbOccCreaMaille
        grCellName = listCreaOccGrName(iOcc)
        nbCell     = listCreaOccNb(iOcc)
        call jeexin(jexnom(meshOut//'.GROUPEMA', grCellName), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(meshOut//'.GROUPEMA', grCellName))
        else
            call utmess('F', 'ALGELINE4_9', sk = grCellName)
        endif
        call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONMAX', max(nbCell , 1))
        call jeecra(jexnom(meshOut//'.GROUPEMA', grCellName), 'LONUTI', nbCell )
        call jeveuo(jexnom(meshOut//'.GROUPEMA', grCellName), 'E', iagma)
        do iCell = 1, nbCell 
            cellName = listCreaName(iCell+shiftCell)
            call jenonu(jexnom(nommai, cellName), zi(iagma-1+iCell))
        end do
        shiftCell = shiftCell + nbCell
    end do
    ASSERT(shiftCell .eq. nbCellCrea)

    call jeexin(grpnov, iret)
    if (iret .eq. 0) then
        nbgrno=0
    else
        call jelira(grpnov, 'NOMUTI', nbgrno)
        call jecreo(gpptnn, 'G N K24')
        call jeecra(gpptnn, 'NOMMAX', nbgrno, ' ')
        call jecrec(grpnoe, 'G V I', 'NO '//gpptnn, 'DISPERSE', 'VARIABLE',&
                    nbgrno)
        do i = 1, nbgrno
            call jenuno(jexnum(grpnov, i), grNodeName)
            call jeveuo(jexnum(grpnov, i), 'L', jvg)
            call jelira(jexnum(grpnov, i), 'LONUTI', nbno)
            call jeexin(jexnom(grpnoe, grCellName), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpnoe, grNodeName))
            else
                call utmess('F', 'ALGELINE4_11', sk=grNodeName)
            endif
            call jeecra(jexnom(grpnoe, grNodeName), 'LONMAX', max(nbno, 1))
            call jeecra(jexnom(grpnoe, grNodeName), 'LONUTI', nbno)
            call jeveuo(jexnom(grpnoe, grNodeName), 'E', jgg)
            do j = 0, nbno-1
                zi(jgg+j)=zi(jvg+j)
            end do
        end do
    endif
!
! - Add new nodes from modified cells
!
    if (nbCellModi .ne. 0) then
        call cmmoma(meshOut, nbCellModi, modiCellNume, modiCellType, nbNodeIn)
    endif
!
! ----------------------------------------------------------------------
!         CREATION DES GROUP_MA ASSOCIE AU MOT CLE "CREA_POI1"
! ----------------------------------------------------------------------
!
    if (nbOccCreaPoi1 .ne. 0) then
        nbOccCreaMaille=0
        do iocc = 1, nbOccCreaPoi1
            call getvtx('CREA_POI1', 'NOM_GROUP_MA', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) nbOccCreaMaille=nbOccCreaMaille+1
        end do
        if (nbOccCreaMaille .ne. 0) then
            call jeexin(grpmai, iret)
            if (iret .eq. 0) then
                call jecreo(gpptnm, 'G N K24')
                call jeecra(gpptnm, 'NOMMAX', nbOccCreaMaille, ' ')
                call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE',&
                            nbOccCreaMaille)
            else
                grpmav='&&OP0167.GROUPEMA'
                call jelira(grpmai, 'NOMUTI', nbgma)
                nbgrmt=nbgma+nbOccCreaMaille
                call cpclma(meshOut, '&&OP0167', 'GROUPEMA', 'V')
                call jedetr(grpmai)
                call jedetr(gpptnm)
                call jecreo(gpptnm, 'G N K24')
                call jeecra(gpptnm, 'NOMMAX', nbgrmt, ' ')
                call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE',&
                            nbgrmt)
                do i = 1, nbgma
                    call jenuno(jexnum(grpmav, i), grCellName)
                    call jeexin(jexnom(grpmai, grCellName), iret)
                    if (iret .eq. 0) then
                        call jecroc(jexnom(grpmai, grCellName))
                    else
                        valk(1)=grCellName
                        call utmess('F', 'ALGELINE4_9', sk=valk(1))
                    endif
                    call jeveuo(jexnum(grpmav, i), 'L', jvg)
                    call jelira(jexnum(grpmav, i), 'LONMAX', nbma)
                    call jeecra(jexnom(grpmai, grCellName), 'LONMAX', max(1, nbma))
                    call jelira(jexnum(grpmav, i), 'LONUTI', nbma)
                    call jeecra(jexnom(grpmai, grCellName), 'LONUTI', nbma)
                    call jeveuo(jexnom(grpmai, grCellName), 'E', jgg)
                    do j = 0, nbma-1
                        zi(jgg+j)=zi(jvg+j)
                    end do
                end do
            endif
            do iocc = 1, nbOccCreaPoi1
                call getvtx('CREA_POI1', 'NOM_GROUP_MA', iocc=iocc, nbval=0, nbret=n1)
                if (n1 .ne. 0) then
                    call getvtx('CREA_POI1', 'NOM_GROUP_MA', iocc=iocc, scal=nogma, nbret=n1)
                    call jenonu(jexnom(grpmai, nogma), ibid)
                    if (ibid .gt. 0) then
                        call utmess('F', 'ALGELINE3_7', sk=nogma)
                    endif
                    call reliem(' ', meshIn, 'NO_NOEUD', motfac, iocc,&
                                nbmc, motcle, tymocl, nomjv, nbma)
                    call jeveuo(nomjv, 'L', jmail)
!
                    call jeexin(jexnom(grpmai, nogma), iret)
                    if (iret .eq. 0) then
                        call jecroc(jexnom(grpmai, nogma))
                    else
                        valk(1)=nogma
                        call utmess('F', 'ALGELINE4_9', sk=valk(1))
                    endif
                    call jeecra(jexnom(grpmai, nogma), 'LONMAX', max(nbma, 1))
                    call jeecra(jexnom(grpmai, nogma), 'LONUTI', nbma)
                    call jeveuo(jexnom(grpmai, nogma), 'E', iagma)
                    do iCell = 0, nbma-1
                        call jenonu(jexnom(nommai, zk8(jmail+iCell)), zi( iagma+iCell))
                    end do
                    if (niv .ge. 1) then
                        write (ifm,902)iocc
                        write (ifm,903)nogma,nbma
                    endif
                endif
            end do
        endif
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
    AS_DEALLOCATE(vk8 = addNodeName)
    AS_DEALLOCATE(vi = addNodeNume)
    AS_DEALLOCATE(vi = listCreaOccNb)
    AS_DEALLOCATE(vk24 = listCreaOccGrName)
!
    call jedema()
!
902 format ('MOT CLE FACTEUR "CREA_POI1", OCCURRENCE ',i4)
903 format ('  CREATION DU GROUP_MA ',a8,' DE ',i6,' MAILLES POI1')
!
end subroutine
