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
use mesh_module, only : checkInclude
!
implicit none
!
#include "asterf_types.h"
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
#include "asterfort/cocali.h"
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
#include "asterfort/lxlgut.h"
#include "asterfort/palim2.h"
#include "asterfort/palim3.h"
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
    integer :: i, lgno, lgnu, nbmc, iret, iad, nbma, iqtr
    integer :: n1, nbrest
    integer :: n1b, jlima, ndinit
    parameter(nbmc=5)
    real(kind=8) :: epais
    character(len=4) :: answer
    character(len=4) :: cdim
    character(len=8) :: prefix
    character(len=8) :: meshIn, meshOut, newmai, model, geofi
    character(len=8) :: nomori, knume, plan, trans
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
    character(len=24) :: momanu, momano, crgrnu, crgrno, lisi
    character(len=24) :: lisk
    character(len=24) :: nomg, valk(2), nogma, gpptnm, gpptnn
    character(len=24) :: prfn1, prfn2, nume2, iadr, nume1, momuto, prfn
    integer :: iaa, iagma, iatyma, ii, ima, in, ino, inumol, j
    integer :: jcrgno, jcrgnu, jgg, jlii, jlik, jmail
    integer :: jmomtu, jnoeu, jnono, jnpt, jopt, jtom, jtrno, jvale, jvg, kvale
    integer :: nbcrp1, nbgma, nbgrma, nbgrmn, nbgrmt, nbgrmv
    integer :: nbgrno, nbmain, nbmaj2, nbmaj3, nbno, nbnot
    integer :: nbpt, nbptt, nori, nrep, ntab, ntpoi
    integer :: ibid, ifm, jdime, jiad, jmomno, jmomnu
    integer :: jnu2, jnum, jpr2, jpro, jrefe, jtypmv
    integer :: nbmaiv, nbmoma, nbnoaj, nbnoev, niv, k, jgeofi
    integer :: dimcon, decala
    integer :: iocc, nbOcc
    character(len=24), parameter :: jvCellNume = '&&OP0167.LISTCELL'
    integer :: nbCell
    integer :: nbField
    integer :: nbOccDecoupeLac, nbOccEclaPg, nbGeomFibre, nbOccCreaFiss, nbOccLineQuad
    integer :: nbOccQuadLine, nbOccModiMaille, nbOccCoquVolu
    integer :: iOccQuadTria
    real(kind=8) :: shrink, lonmin
    aster_logical :: lpb, l_modi_maille
    integer :: prefNume
    character(len=8) :: prefCell, prefNode
    integer, pointer :: listCellNume(:) => null()
    character(len=16), pointer :: listField(:) => null()
    integer, pointer :: adrjvx(:) => null()
    integer, pointer :: nbnoma(:) => null()
    integer, pointer :: nbnomb(:) => null()
    integer, pointer :: nomnum(:) => null()
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
    call getfac('DECOUPE_LAC', nbOccDecoupeLac)
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
        call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNode)
        call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNume)
        call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
        if (nbCell .ne. nbmaiv) then
            call utmess('A', 'MESH1_4', sk=keywfact)
        endif
        call jeveuo(jvCellNume, 'L', vi = listCellNume)
        call cmlqlq(meshIn, meshOut, nbCell, listCellNume, prefNode, prefNume)
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
            call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNode)
            call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNume)
!
            call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
            call jeveuo(jvCellNume, 'L', vi = listCellNume)
            if (nbCell .ne. nbmaiv) then
                call utmess('A', 'MESH1_4', sk=keywfact)
            endif
!
            if (keywfact .eq. 'HEXA20_27') then
                call cm2027(meshIn, meshOut, nbCell, listCellNume, prefNode, prefNume)
            else if (keywfact.eq.'PENTA15_18') then
                call cm1518(meshIn, meshOut, nbCell, listCellNume, prefNode, prefNume)
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
        call getvtx("MODI_HHO", 'PREF_NOEUD', iocc=1, scal=prefix, nbret=n1)
        call getvis("MODI_HHO", 'PREF_NUME', iocc=1, scal=ndinit, nbret=n1)
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
        call cmhho(meshIn, meshOut, nbma, zi(jlima), prefix, ndinit)
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
            call getvtx(keywfact, 'PREF_MAILLE', iocc=iOccQuadTria, scal=prefCell)
            call getvis(keywfact, 'PREF_NUME', iocc=iOccQuadTria, scal=prefNume)
            call getelem(meshIn, keywfact, iOccQuadTria, 'F', jvCellNume, nbCell)
            if (nbCell .ne. nbmaiv) then
                call utmess('A', 'MESH1_4', sk=keywfact)
            endif
            call jeveuo(jvCellNume, 'L', vi = listCellNume)
            call cmqutr('G', meshIn, meshOut, nbCell, listCellNume,&
                        prefCell, prefNume)
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
        call getvtx(keywfact, 'PREF_NOEUD', iocc=1, scal=prefNode)
        call getvtx(keywfact, 'PREF_MAILLE', iocc=1, scal=prefCell)
        call getvis(keywfact, 'PREF_NUME', iocc=1, scal=prefNume)
        call getvtx(keywfact, 'PLAN', iocc=1, scal=plan)
        if (plan .eq. 'MOY') then
            trans = 'INF'
            call getvtx(keywfact, 'TRANSLATION', iocc=1, scal=trans)
        endif
        call getelem(meshIn, keywfact, 1, 'F', jvCellNume, nbCell)
        call jeveuo(jvCellNume, 'L', vi = listCellNume)
        call cmcovo(meshIn, meshOut, nbCell, listCellNume, prefNode,&
                    prefCell, prefNume, epais, plan, trans)
        goto 350
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "RESTREINT"
!
! --------------------------------------------------------------------------------------------------
!
    call getfac('RESTREINT', nbrest)
    if (nbrest .ne. 0) then
        call rdtmai(meshIn, meshOut, 'G', meshOut//'.CRNO', meshOut// '.CRMA',&
                    'G', 0, [0])
! ---    VERIFICATIONS DU MAILLAGE
        call chckma(meshOut, 1.0d-03)
        goto 350
!
    endif
!
! ----------------------------------------------------------------------
!               AURES MOTS CLES :
! ----------------------------------------------------------------------
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
    nbnoev=zi(jdime)
    nbmaiv=zi(jdime+3-1)
!
    call jeveuo(typmav, 'L', jtypmv)
!
! --------------------------------------------------------------------------------------------------
!
!   For "MODI_MAILLE"
!
! --------------------------------------------------------------------------------------------------
!
    call getfac('MODI_MAILLE', nbmoma)
    nbnoaj=0
    if (nbmoma .ne. 0) then
        momanu='&&OP0167.MO_MA.NUM'
        momano='&&OP0167.MO_MA.NOM'
!
        momuto='&&OP0167.MO_TO.NUM'
!
        lisi='&&OP0167.LISI'
        lisk='&&OP0167.LISK'
!
        iadr='&&OP0167.IADR'
        prfn='&&OP0167.PRFN'
        nume1='&&OP0167.NUME'
        prfn2='&&OP0167.PRFN2'
        nume2='&&OP0167.NUME2'
!
        call wkvect(momanu, 'V V I', nbmaiv, jmomnu)
        call wkvect(momano, 'V V K8', nbmaiv, jmomno)
!
        call wkvect(iadr, 'V V I', nbmoma, jiad)
        call wkvect(prfn, 'V V K8', nbmoma, jpro)
        call wkvect(nume1, 'V V I', nbmoma, jnum)
        call wkvect(prfn2, 'V V K8', nbmaiv, jpr2)
        call wkvect(nume2, 'V V I', nbmaiv, jnu2)
!
        l_modi_maille = ASTER_FALSE
!
        iad=1
        do iocc = 1, nbmoma
            call getvtx('MODI_MAILLE', 'OPTION', iocc=iocc, scal=option, nbret=n1)
            zi(jiad+iocc-1)=1
            call getvtx('MODI_MAILLE', 'PREF_NOEUD', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvtx('MODI_MAILLE', 'PREF_NOEUD', iocc=iocc, scal=zk8(jpro+iocc-1),&
                            nbret=n1)
                lgno=lxlgut(zk8(jpro+iocc-1))
            endif
            call getvis('MODI_MAILLE', 'PREF_NUME', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvis('MODI_MAILLE', 'PREF_NUME', iocc=iocc, scal=zi(jnum+iocc-1),&
                            nbret=n1)
            endif
            call palim2('MODI_MAILLE', iocc, meshIn, momanu, momano,&
                        zi(jiad+iocc-1))
            if (zi(jiad+iocc-1)-1 .le. 0) then
                call utmess('A', 'MODELISA3_32', sk=option, si=iocc)
                goto 60
            else
                l_modi_maille = ASTER_TRUE
            endif
!
            call wkvect(lisi, 'V V I', zi(jiad+iocc-1)-1, jlii)
            call wkvect(lisk, 'V V K8', zi(jiad+iocc-1)-1, jlik)
!
            do ii = 1, zi(jiad+iocc-1)-1
                zi(jlii+ii-1)=zi(jmomnu+ii-1)
                zk8(jlik+ii-1)=zk8(jmomno+ii-1)
            end do
            call cocali(momuto, lisi, 'I')
            iaa=iad
            iad=iad+zi(jiad+iocc-1)-1
!
! LE PREFIXE EST LE MEME POUR TOUS LES NOEUDS ENTRE
! L'ANCIENNE ET LA NOUVELLE ADRESSE
!
            do ii = iaa, iad-1
                zk8(jpr2+ii-1)=zk8(jpro+iocc-1)
            end do
!
! LE PREF_NUME EST A DEFINIR POUR LE PREMIER NOEUD
! LES AUTRES SE TROUVENT EN INCREMENTANT
!
            zi(jnu2+iaa-1)=zi(jnum+iocc-1)
            call jedetr(lisi)
            call jedetr(lisk)
!
            if (niv .ge. 1) then
                write (ifm,900)iocc
                if (option .eq. 'TRIA6_7') then
                    write (ifm,901)zi(jiad+iocc-1)-1,'TRIA6','TRIA7'
                else if (option.eq.'QUAD8_9') then
                    write (ifm,901)zi(jiad+iocc-1)-1,'QUAD8','QUAD9'
                else if (option.eq.'SEG3_4') then
                    write (ifm,901)zi(jiad+iocc-1)-1,'SEG3','SEG4'
                endif
            endif
 60         continue
        end do
!
        if(l_modi_maille) then
            call jeveuo(momuto, 'L', jmomtu)
        else
            call utmess('A', 'MODELISA3_31', si=nbmoma)
            nbmoma = 0
        end if
!
        nbnoaj=iad-1
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "CREA_MAILLE"
!
! --------------------------------------------------------------------------------------------------
!
    call getfac('CREA_MAILLE', nbgrma)
    nbmaj2=0
    if (nbgrma .ne. 0) then
        crgrnu='&&OP0167.CR_GR.NUM'
        crgrno='&&OP0167.CR_GR.NOM'
        call wkvect(crgrnu, 'V V I', nbmaiv, jcrgnu)
        call wkvect(crgrno, 'V V K8', nbmaiv, jcrgno)
        nbmaj2=0
        do iocc = 1, nbgrma
            call palim3('CREA_MAILLE', iocc, meshIn, crgrnu, crgrno,&
                        nbmaj2)
        end do
        call jeveuo(crgrnu, 'L', jcrgnu)
        call jeveuo(crgrno, 'L', jcrgno)
    endif
!
! --------------------------------------------------------------------------------------------------
!
!   For "CREA_POI1"
!
! --------------------------------------------------------------------------------------------------
!
    call getfac('CREA_POI1', nbcrp1)
    nbmaj3=0
    if (nbcrp1 .ne. 0) then
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
        call wkvect('&&OP0167.IND_NOEUD', 'V V I', nbnoev, jtrno)
        call wkvect('&&OP0167.NOM_NOEUD', 'V V K8', nbnoev, jnono)
!
        do iocc = 1, nbcrp1
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
        do ima = 1, nbnoev
            if (zi(jtrno+ima-1) .eq. 0) goto 110
            call jenuno(jexnum(nomnov, ima), newmai)
            call jenonu(jexnom(nommav, newmai), ibid)
            if (ibid .eq. 0) then
                nbmaj3=nbmaj3+1
                zk8(jnono-1+nbmaj3)=newmai
            else
                valk(1)=newmai
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
    if (nbnoaj .ne. 0) then
        nbnot=nbnoev+nbnoaj
        zi(jdime)=nbnot
!
        call jecreo(nomnoe, 'G N K8')
        call jeecra(nomnoe, 'NOMMAX', nbnot, ' ')
        do ino = 1, nbnoev
            call jenuno(jexnum(nomnov, ino), nomg)
            call jeexin(jexnom(nomnoe, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(nomnoe, nomg))
            else
                valk(1)=nomg
                call utmess('F', 'ALGELINE4_5', sk=valk(1))
            endif
        end do
        do ino = nbnoev+1, nbnot
! TRAITEMENT DES NOEUDS AJOUTES
! ON CODE LE NUMERO DU NOEUD COURANT
            call codent(zi(jnu2+ino-nbnoev-1), 'G', knume)
!
! SI LE PREFIXE COURANT EST LE MEME QUE LE SUIVANT ALORS
! LE NUME EST INCREMENTE
            if (zk8(jpr2+ino-nbnoev-1) .eq. zk8(jpr2+ino-nbnoev)) then
                zi(jnu2+ino-nbnoev)=zi(jnu2+ino-nbnoev-1)+1
            endif
!
            lgnu=lxlgut(knume)
            prfn1=zk8(jpr2+ino-nbnoev-1)
            lgno=lxlgut(prfn1)
            if (lgnu+lgno .gt. 8) then
                call utmess('F', 'ALGELINE_16')
            endif
            nomg=prfn1(1:lgno)//knume
            call jeexin(jexnom(nomnoe, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(nomnoe, nomg))
            else
                valk(1)=nomg
                call utmess('F', 'ALGELINE4_5', sk=valk(1))
            endif
        end do
!
        call jeveuo(coovav, 'L', jvale)
        call wkvect(cooval, 'G V R8', 3*nbnot, kvale)
        do i = 0, 3*nbnoev-1
            zr(kvale+i)=zr(jvale+i)
        end do
        call jelira(coovav, 'DOCU', cval=cdim)
        call jeecra(cooval, 'DOCU', cval=cdim)
    else
        call jedupo(nomnov, 'G', nomnoe, .false._1)
        call jedupo(coovav, 'G', cooval, .false._1)
    endif
!
! --- CAS OU L'ON FOURNIT UNE TABLE.
! --- IL S'AGIT DE DEFINIR LES COORDONNEES DES NOEUDS DU MAILLAGE
! --- EN SORTIE DANS UN NOUVEAU REPERE.
! --- CETTE FONCTIONNALITE SERT DANS LE CAS OU L'ON CALCULE LES
! --- CARACTERISTIQUES DE CISAILLEMENT D'UNE POUTRE A PARTIR DE LA
! --- DONNEE D'UNE SECTION DE CETTE POUTRE MAILLEE AVEC DES ELEMENTS
! --- MASSIFS 2D.
! --- LA TABLE OBTENUE PAR POST_ELEM (OPTION : CARA_GEOM)  CONTIENT
! --- LES COORDONNEES DE LA NOUVELLE ORIGINE  (I.E. LE CENTRE DE
! --- GRAVITE) ET L'ANGLE FORME PAR LES AXES PRINCIPAUX D'INERTIE
! --- (LES NOUVEAUX AXES) AVEC LES AXES GLOBAUX :
! --- ON DEFINIT LE MAILLAGE EN SORTIE DANS CE NOUVEAU REPERE
! --- POUR LE CALCUL DU CENTRE DE CISAILLEMENT TORSION ET DES
! --- COEFFICIENTS DE CISAILLEMENT.
! --- DANS LE CAS OU L'ON DONNE LE MOT-CLE ORIG_TORSION
! --- LA TABLE CONTIENT LES COORDONNEES DU CENTRE DE CISAILLEMENT-
! --- TORSION ET ON DEFINIT LE NOUVEAU MAILLAGE EN PRENANT COMME
! --- ORIGINE CE POINT. CETTE OPTION EST UTILISEE POUR LE CALCUL
! --- DE L'INERTIE DE GAUCHISSEMENT :
!     -----------------------------
    call getfac('REPERE', nrep)
    if (nrep .ne. 0) then
        call getvid('REPERE', 'TABLE', iocc=1, nbval=0, nbret=ntab)
        if (ntab .ne. 0) then
            call getvid('REPERE', 'TABLE', iocc=1, scal=table, nbret=ntab)
            call getvtx('REPERE', 'NOM_ORIG', iocc=1, nbval=0, nbret=nori)
            if (nori .ne. 0) then
                call getvtx('REPERE', 'NOM_ORIG', iocc=1, scal=nomori, nbret=nori)
                if (nomori .eq. 'CDG') then
                    call chcoma(table, meshOut)
                else if (nomori.eq.'TORSION') then
                    call chcomb(table, meshOut)
                else
                    call utmess('F', 'ALGELINE3_5')
                endif
            endif
        endif
    endif
!
! ----------------------------------------------------------------------
!         ON AGRANDIT LE '.NOMMAI' ET LE '.CONNEX'
! ----------------------------------------------------------------------
!
    nbmain=nbmaiv+nbmaj2+nbmaj3
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
    do ima = 1, nbmaiv
        call jenuno(jexnum(nommav, ima), nomg)
        call jeexin(jexnom(nommai, nomg), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(nommai, nomg))
        else
            valk(1)=nomg
            call utmess('F', 'ALGELINE4_7', sk=valk(1))
        endif
!
        call jenonu(jexnom(nommav, nomg), ibid)
        jtom=jtypmv-1+ibid
        call jenonu(jexnom(nommai, nomg), ibid)
        zi(iatyma-1+ibid)=zi(jtom)
!
        call jenonu(jexnom(nommav, nomg), ibid)
        call jelira(jexnum(connev, ibid), 'LONMAX', nbpt)
        call jeveuo(jexnum(connev, ibid), 'L', jopt)
        nbptt=nbpt
        do in = 1, nbnoaj
            if (ima .eq. zi(jmomtu+in-1)) then
                nbptt=nbpt+1
                goto 160
!
            endif
        end do
160     continue
        call jenonu(jexnom(nommai, nomg), ibid)
        dimcon = dimcon+nbptt
        nbnoma(ima) = nbptt
        nbnomb(ima) = nbpt
        adrjvx(ima) = jopt
        nomnum(ima) = ibid
    end do
!
    decala = decala + nbmaiv
!
    do ima = 1, nbmaj2
        newmai=zk8(jcrgno+ima-1)
        inumol=zi(jcrgnu+ima-1)
        call jeexin(jexnom(nommai, newmai), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(nommai, newmai))
        else
            valk(1)=newmai
            call utmess('F', 'ALGELINE4_7', sk=valk(1))
        endif
!
        jtom=jtypmv-1+inumol
        call jenonu(jexnom(nommai, newmai), ibid)
        if (ibid .eq. 0) then
            call utmess('F', 'ALGELINE3_6', sk=newmai)
        endif
        zi(iatyma-1+ibid)=zi(jtom)
!
        call jelira(jexnum(connev, inumol), 'LONMAX', nbpt)
        call jeveuo(jexnum(connev, inumol), 'L', jopt)
        dimcon = dimcon+nbpt
        nbnoma(1+decala+ima-1) = nbpt
        adrjvx(1+decala+ima-1) = jopt
        nomnum(1+decala+ima-1) = ibid
    end do
!
    dimcon = dimcon+nbmaj3
    call jeecra(connex, 'LONT', dimcon)
!
    decala = 0
    do ima = 1, nbmaiv
        nbptt = nbnoma(1+decala+ima-1)
        nbpt = nbnomb(1+decala+ima-1)
        jopt = adrjvx(1+decala+ima-1)
        ibid = nomnum(1+decala+ima-1)
        call jeecra(jexnum(connex, ibid), 'LONMAX', nbptt)
        call jeveuo(jexnum(connex, ibid), 'E', jnpt)
        do ino = 0, nbpt-1
            zi(jnpt+ino)=zi(jopt+ino)
        end do
    end do
!
    decala = decala + nbmaiv
!
    do ima = 1, nbmaj2
        nbpt = nbnoma(1+decala+ima-1)
        jopt = adrjvx(1+decala+ima-1)
        ibid = nomnum(1+decala+ima-1)
        call jeecra(jexnum(connex, ibid), 'LONMAX', nbpt)
        call jeveuo(jexnum(connex, ibid), 'E', jnpt)
        do ino = 0, nbpt-1
            zi(jnpt+ino)=zi(jopt+ino)
        end do
    end do
!
    do ima = 1, nbmaj3
        newmai=zk8(jnono+ima-1)
        call jenonu(jexnom(nommai, newmai), ibid)
        if (ibid .ne. 0) goto 230
        call jeexin(jexnom(nommai, newmai), iret)
        if (iret .eq. 0) then
            call jecroc(jexnom(nommai, newmai))
        else
            valk(1)=newmai
            call utmess('F', 'ALGELINE4_7', sk=valk(1))
        endif
!
        call jenonu(jexnom(nommai, newmai), ibid)
        if (ibid .eq. 0) then
            call utmess('F', 'ALGELINE3_6', sk=newmai)
        endif
        zi(iatyma-1+ibid)=ntpoi
!
        call jeecra(jexnum(connex, ibid), 'LONMAX', 1)
        call jeveuo(jexnum(connex, ibid), 'E', jnpt)
        call jenonu(jexnom(nomnoe, newmai), zi(jnpt))
230     continue
    end do
    AS_DEALLOCATE(vi=nbnoma)
    AS_DEALLOCATE(vi=nbnomb)
    AS_DEALLOCATE(vi=adrjvx)
    AS_DEALLOCATE(vi=nomnum)
! ----------------------------------------------------------------------
!
    call jeexin(grpmav, iret)
    if (iret .eq. 0) then
        nbgrmv=0
    else
        call jelira(grpmav, 'NOMUTI', nbgrmv)
    endif
    nbgrmn=nbgrmv+nbgrma
    if (nbgrmn .ne. 0) then
        call jecreo(gpptnm, 'G N K24')
        call jeecra(gpptnm, 'NOMMAX', nbgrmn, ' ')
        call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE',&
                    nbgrmn)
        do i = 1, nbgrmv
            call jenuno(jexnum(grpmav, i), nomg)
            call jeexin(jexnom(grpmai, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpmai, nomg))
            else
                valk(1)=nomg
                call utmess('F', 'ALGELINE4_9', sk=valk(1))
            endif
            call jeveuo(jexnum(grpmav, i), 'L', jvg)
            call jelira(jexnum(grpmav, i), 'LONMAX', nbma)
            call jeecra(jexnom(grpmai, nomg), 'LONMAX', max(nbma, 1))
            call jelira(jexnum(grpmav, i), 'LONUTI', nbma)
            call jeecra(jexnom(grpmai, nomg), 'LONUTI', nbma)
            call jeveuo(jexnom(grpmai, nomg), 'E', jgg)
            do j = 0, nbma-1
                zi(jgg+j)=zi(jvg+j)
            end do
        end do
        do i = 1, nbgrma
            call getvtx('CREA_MAILLE', 'NOM', iocc=i, scal=nomg, nbret=n1)
            ASSERT(n1.eq.1)
            call jeexin(jexnom(grpmai, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpmai, nomg))
            else
                valk(1)=nomg
                call utmess('F', 'ALGELINE4_9', sk=valk(1))
            endif
            nbmaj2=0
            call palim3('CREA_MAILLE', i, meshIn, crgrnu, crgrno,&
                        nbmaj2)
            call jeveuo(crgrno, 'L', jcrgno)
            call jeecra(jexnom(grpmai, nomg), 'LONMAX', max(nbmaj2, 1))
            call jeecra(jexnom(grpmai, nomg), 'LONUTI', nbmaj2)
            call jeveuo(jexnom(grpmai, nomg), 'E', iagma)
            do ima = 0, nbmaj2-1
                call jenonu(jexnom(nommai, zk8(jcrgno+ima)), zi(iagma+ ima))
            end do
        end do
    endif
!
! ----------------------------------------------------------------------
!
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
            call jenuno(jexnum(grpnov, i), nomg)
            call jeveuo(jexnum(grpnov, i), 'L', jvg)
            call jelira(jexnum(grpnov, i), 'LONUTI', nbno)
            call jeexin(jexnom(grpnoe, nomg), iret)
            if (iret .eq. 0) then
                call jecroc(jexnom(grpnoe, nomg))
            else
                valk(1)=nomg
                call utmess('F', 'ALGELINE4_11', sk=valk(1))
            endif
            call jeecra(jexnom(grpnoe, nomg), 'LONMAX', max(nbno, 1))
            call jeecra(jexnom(grpnoe, nomg), 'LONUTI', nbno)
            call jeveuo(jexnom(grpnoe, nomg), 'E', jgg)
            do j = 0, nbno-1
                zi(jgg+j)=zi(jvg+j)
            end do
        end do
    endif
!
    if (nbmoma .ne. 0) call cmmoma(meshOut, momuto, nbnoev, nbnoaj)
!
!
! ----------------------------------------------------------------------
!         CREATION DES GROUP_MA ASSOCIE AU MOT CLE "CREA_POI1"
! ----------------------------------------------------------------------
!
    if (nbcrp1 .ne. 0) then
        nbgrma=0
        do iocc = 1, nbcrp1
            call getvtx('CREA_POI1', 'NOM_GROUP_MA', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) nbgrma=nbgrma+1
        end do
        if (nbgrma .ne. 0) then
            call jeexin(grpmai, iret)
            if (iret .eq. 0) then
                call jecreo(gpptnm, 'G N K24')
                call jeecra(gpptnm, 'NOMMAX', nbgrma, ' ')
                call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE',&
                            nbgrma)
            else
                grpmav='&&OP0167.GROUPEMA'
                call jelira(grpmai, 'NOMUTI', nbgma)
                nbgrmt=nbgma+nbgrma
                call cpclma(meshOut, '&&OP0167', 'GROUPEMA', 'V')
                call jedetr(grpmai)
                call jedetr(gpptnm)
                call jecreo(gpptnm, 'G N K24')
                call jeecra(gpptnm, 'NOMMAX', nbgrmt, ' ')
                call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE',&
                            nbgrmt)
                do i = 1, nbgma
                    call jenuno(jexnum(grpmav, i), nomg)
                    call jeexin(jexnom(grpmai, nomg), iret)
                    if (iret .eq. 0) then
                        call jecroc(jexnom(grpmai, nomg))
                    else
                        valk(1)=nomg
                        call utmess('F', 'ALGELINE4_9', sk=valk(1))
                    endif
                    call jeveuo(jexnum(grpmav, i), 'L', jvg)
                    call jelira(jexnum(grpmav, i), 'LONMAX', nbma)
                    call jeecra(jexnom(grpmai, nomg), 'LONMAX', max(1, nbma))
                    call jelira(jexnum(grpmav, i), 'LONUTI', nbma)
                    call jeecra(jexnom(grpmai, nomg), 'LONUTI', nbma)
                    call jeveuo(jexnom(grpmai, nomg), 'E', jgg)
                    do j = 0, nbma-1
                        zi(jgg+j)=zi(jvg+j)
                    end do
                end do
            endif
            do iocc = 1, nbcrp1
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
                    do ima = 0, nbma-1
                        call jenonu(jexnom(nommai, zk8(jmail+ima)), zi( iagma+ima))
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
350 continue

! - Add title in mesh datastructure
    call titre()

! - Update parameters for modified mesh (bounding box and dimensions)
    call cargeo(meshOut)

! - Verbose
    call infoma(meshOut)
!
    call jedema()
!
900 format ('MOT CLE FACTEUR "MODI_MAILLE", OCCURRENCE ',i4)
901 format ('  MODIFICATION DE ',i6,' MAILLES ',a8,' EN ',a8)
902 format ('MOT CLE FACTEUR "CREA_POI1", OCCURRENCE ',i4)
903 format ('  CREATION DU GROUP_MA ',a8,' DE ',i6,' MAILLES POI1')
!
end subroutine
