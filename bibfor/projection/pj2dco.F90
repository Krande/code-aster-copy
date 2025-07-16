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
subroutine pj2dco(typeSelect, &
                  entity1, entity2, &
                  nbCellSelect1, listCellSelect1, &
                  nbNodeSelect2, listNodeSelect2, &
                  geom1, geom2, corrMesh, &
                  l_dmax, dmax, dala, &
                  listInterc_, nbInterc_)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/pj2dap.h"
#include "asterfort/pj2dfb.h"
#include "asterfort/pj2dtr.h"
#include "asterfort/pjxxut.h"
#include "asterfort/utimsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*), intent(in) :: typeSelect
    character(len=8), intent(in) :: entity1, entity2
    integer(kind=8), intent(in) :: nbCellSelect1, listCellSelect1(*)
    integer(kind=8), intent(in) :: nbNodeSelect2, listNodeSelect2(*)
    character(len=*), intent(in) :: geom1, geom2
    character(len=16), intent(in)  :: corrMesh
    aster_logical, intent(in) :: l_dmax
    real(kind=8), intent(in) :: dmax, dala
    character(len=16), optional, intent(in)  :: listInterc_
    integer(kind=8), optional, intent(in)  :: nbInterc_
!
! --------------------------------------------------------------------------------------------------
!
! Create corresp_2_mailla datastructure
!
! Projection of nodes from entity1 to entity2 when entity1 is 2D
!
! --------------------------------------------------------------------------------------------------
!
! In  typeSelect       : type of selection (all or restricted) TOUT or PARTIE
! In  entity1          : name of first entity (model or mesh)
! In  entity2          : name of second entity (model or mesh)
! In  nbCellSelect1    : number of cells to select in first entity if typeSelect = 'PARTIE'
! In  listCellSelect1  : list of cells to select in first entity if typeSelect = 'PARTIE'
! In  nbNodeSelect2    : number of nodes to select in second entity if typeSelect = 'PARTIE'
! In  listNodeSelect2  : list of nodes to select in second entity if typeSelect = 'PARTIE'
! In  geom1            : name of object with coordinates of nodes for entity 1
! In  geom2            : name of object with coordinates of nodes for entity 2
! In  coorMesh         : name of corresp_2_mailla datastructure
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: dbg = ASTER_FALSE
    integer(kind=8), parameter :: spacedim = 2
    character(len=8) :: mesh1, mesh2, nodeName2
    character(len=14), parameter :: boite = '&&PJ2DCO.BOITE'
    character(len=16), parameter :: corrMeshTemp = '&&PJ2DCO.CORRESP'
    character(len=16) :: listInterc
    integer(kind=8) :: nbCellType
    integer(kind=8) :: cellListType(MT_NTYMAX)
    character(len=8) :: cellListCode(MT_NTYMAX)
    integer(kind=8) :: ifm, niv, nbNode1, nbNode2, nbCell1, nbCell2, k
    integer(kind=8) :: nbTria, iatr3
    integer(kind=8) :: iCell1, iNode2
    integer(kind=8) :: iacoo1, iacoo2, nbpt0, ino2_0, idecal_0
    integer(kind=8) :: iabtco, jxxk1, iaconu, iacocf, iacotr
    integer(kind=8) :: ilcnx1
    integer(kind=8) :: iaconb, cellTypeNume, idecal, cellLink1, nbtrou, nbInterc
    aster_logical :: loin
    real(kind=8) :: dmin, cobary(3)
    integer(kind=8), pointer :: bt2dlc(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: bt2ddi(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: bt2dnb(:) => null()
    real(kind=8), pointer :: bt2dvr(:) => null()
    integer(kind=8), pointer :: listCell1(:) => null(), listNode2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

! - Prepare list of entities
    call pjxxut('2D', typeSelect, &
                entity1, entity2, &
                nbCellSelect1, listCellSelect1, &
                nbNodeSelect2, listNodeSelect2, &
                mesh1, mesh2, &
                nbCellType, cellListType, cellListCode)
    call jeveuo('&&PJXXCO.LIMA1', 'L', vi=listCell1)
    call jeveuo('&&PJXXCO.LINO2', 'L', vi=listNode2)

! - Access to meshes
    call dismoi('NB_NO_MAILLA', mesh1, 'MAILLAGE', repi=nbNode1)
    call dismoi('NB_NO_MAILLA', mesh2, 'MAILLAGE', repi=nbNode2)
    call dismoi('NB_MA_MAILLA', mesh1, 'MAILLAGE', repi=nbCell1)
    call dismoi('NB_MA_MAILLA', mesh2, 'MAILLAGE', repi=nbCell2)
    call jeveuo(mesh1//'.TYPMAIL', 'L', vi=typmail)

! - Cut cells in TRIA3: number of TRIA3
    nbTria = 0
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            cellTypeNume = typmail(iCell1)
            if (cellTypeNume .eq. cellListType(1)) then
                nbTria = nbTria+1
            else if (cellTypeNume .eq. cellListType(2)) then
                nbTria = nbTria+1
            else if (cellTypeNume .eq. cellListType(3)) then
                nbTria = nbTria+1
            else if (cellTypeNume .eq. cellListType(4)) then
                nbTria = nbTria+2
            else if (cellTypeNume .eq. cellListType(5)) then
                nbTria = nbTria+2
            else if (cellTypeNume .eq. cellListType(6)) then
                nbTria = nbTria+2
            else
                call utmess('F', 'PROJECTION4_1')
            end if
        end if
    end do
    if (nbTria .eq. 0) then
        call utmess('F', 'PROJECTION4_55')
    end if

! - Create object to cut cells in TRIA3
!         long(v)=1+4*ntr3
!         v(1) : ntr3(=nombre de tria3)
!         v(1+4(i-1)+1) : numero du 1er  noeud du ieme tria3
!         v(1+4(i-1)+2) : numero du 2eme noeud du ieme tria3
!         v(1+4(i-1)+3) : numero du 3eme noeud du ieme tria3
!         v(1+4(i-1)+4) : numero de la maille mere du ieme tria3
!
    call wkvect('&&PJXXCO.TRIA3', 'V V I', 1+4*nbTria, iatr3)
    zi(iatr3-1+1) = nbTria
    call jeveuo(mesh1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    nbTria = 0
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            cellTypeNume = typmail(iCell1)
            if ((cellTypeNume .eq. cellListType(1)) .or. &
                (cellTypeNume .eq. cellListType(2)) .or. &
                (cellTypeNume .eq. cellListType(3))) then
                nbTria = nbTria+1
                zi(iatr3+(nbTria-1)*4+4) = iCell1
                zi(iatr3+(nbTria-1)*4+1) = connex(1+zi(ilcnx1-1+iCell1)-2+1)
                zi(iatr3+(nbTria-1)*4+2) = connex(1+zi(ilcnx1-1+iCell1)-2+2)
                zi(iatr3+(nbTria-1)*4+3) = connex(1+zi(ilcnx1-1+iCell1)-2+3)
            else if ((cellTypeNume .eq. cellListType(4)) .or. &
                     (cellTypeNume .eq. cellListType(5)) .or. &
                     (cellTypeNume .eq. cellListType(6))) then
                nbTria = nbTria+1
                zi(iatr3+(nbTria-1)*4+4) = iCell1
                zi(iatr3+(nbTria-1)*4+1) = connex(1+zi(ilcnx1-1+iCell1)-2+1)
                zi(iatr3+(nbTria-1)*4+2) = connex(1+zi(ilcnx1-1+iCell1)-2+2)
                zi(iatr3+(nbTria-1)*4+3) = connex(1+zi(ilcnx1-1+iCell1)-2+3)
                nbTria = nbTria+1
                zi(iatr3+(nbTria-1)*4+4) = iCell1
                zi(iatr3+(nbTria-1)*4+1) = connex(1+zi(ilcnx1-1+iCell1)-2+1)
                zi(iatr3+(nbTria-1)*4+2) = connex(1+zi(ilcnx1-1+iCell1)-2+3)
                zi(iatr3+(nbTria-1)*4+3) = connex(1+zi(ilcnx1-1+iCell1)-2+4)
            else
                call utmess('F', 'PROJECTION4_1')
            end if
        end if
    end do

! - Get coordinates of nodes
    if (geom1 .eq. ' ') then
        call jeveuo(mesh1//'.COORDO    .VALE', 'L', iacoo1)
    else
        call jeveuo(geom1, 'L', iacoo1)
    end if
    if (geom2 .eq. ' ') then
        call jeveuo(mesh2//'.COORDO    .VALE', 'L', iacoo2)
    else
        call jeveuo(geom2, 'L', iacoo2)
    end if

! - Create boxes
    call pj2dfb(boite, zi(iatr3), zr(iacoo1), zr(iacoo2))
    call jeveuo(boite//'.BT2DDI', 'L', vi=bt2ddi)
    call jeveuo(boite//'.BT2DVR', 'L', vr=bt2dvr)
    call jeveuo(boite//'.BT2DNB', 'L', vi=bt2dnb)
    call jeveuo(boite//'.BT2DLC', 'L', vi=bt2dlc)
    call jeveuo(boite//'.BT2DCO', 'L', iabtco)

! Boxes for 2D case
!     boite_2d (k14) ::= record
!      .bt2ddi   : ojb s v i  long=2
!      .bt2dvr   : ojb s v r  long=6
!      .bt2dnb   : ojb s v i  long=nx*ny
!      .bt2dlc   : ojb s v i  long=1+nx*ny
!      .bt2dco   : ojb s v i  long=*
!      .bt2ddi(1) : nx=nombre de boites dans la direction x
!      .bt2ddi(2) : ny=nombre de boites dans la direction y
!      .bt2dvr(1) : xmin     .bt2dvr(2) : xmax
!      .bt2dvr(3) : ymin     .bt2dvr(4) : ymax
!      .bt2dvr(5) : dx = (xmax-xmin)/nbx
!      .bt2dvr(6) : dy = (ymax-ymin)/nby
!      .bt2dnb    : longueurs des boites
!      .bt2dnb(1) : nombre de tria3 contenus dans la boite(1,1)
!      .bt2dnb(2) : nombre de tria3 contenus dans la boite(2,1)
!      .bt2dnb(3) : ...
!      .bt2dnb(nx*ny) : nombre de tria3 contenus dans la boite(nx,ny)
!      .bt2dlc    : longueurs cumulees de .bt2dco
!      .bt2dlc(1) : 0
!      .bt2dlc(2) : bt2dlc(1)+nbtr3(boite(1,1))
!      .bt2dlc(3) : bt2dlc(2)+nbtr3(boite(2,1))
!      .bt2dlc(4) : ...
!      .bt2dco    : contenu des boites
!       soit   nbtr3 =nbtr3(boite(p,q)=bt2dnb((q-1)*nx+p)
!              debtr3=bt2dlc((q-1)*nx+p)
!        do k=1,nbtr3
!          tr3=.bt2dco(debtr3+k)
!        done
!        tr3 est le numero du kieme tria3 de la boite (p,q)

! - Create temporary datastructure
    call wkvect(corrMeshTemp//'.PJXX_K1', 'V V K24', 5, jxxk1)
    zk24(jxxk1-1+1) = mesh1
    zk24(jxxk1-1+2) = mesh2
    zk24(jxxk1-1+3) = 'COLLOCATION'
    call wkvect(corrMeshTemp//'.PJEF_NB', 'V V I', nbNode2, iaconb)
    call wkvect(corrMeshTemp//'.PJEF_NU', 'V V I', 3*nbNode2, iaconu)
    call wkvect(corrMeshTemp//'.PJEF_CF', 'V V R', 3*nbNode2, iacocf)
    call wkvect(corrMeshTemp//'.PJEF_TR', 'V V I', nbNode2, iacotr)

! - Look for cell cellLink1 in entity1 for each node from entity 2
    idecal = 0
    nbpt0 = 0
    do iNode2 = 1, nbNode2
        if (listNode2(iNode2) .eq. 0) cycle
! ----- Special for XFEM (issue23983))
        if (zr(iacoo2-1+3*(iNode2-1)+1) .eq. 0.d0 .and. &
            zr(iacoo2-1+3*(iNode2-1)+2) .eq. 0.d0) then
            nbpt0 = nbpt0+1
            if (nbpt0 .eq. 1) then
                ino2_0 = iNode2
                idecal_0 = idecal
            else
                zi(iaconb-1+iNode2) = 3
                cellLink1 = zi(iacotr-1+ino2_0)
                zi(iacotr-1+iNode2) = cellLink1
                if (cellLink1 .eq. 0) cycle
                do k = 1, 3
                    zi(iaconu-1+idecal+k) = zi(iaconu-1+idecal_0+k)
                    zr(iacocf-1+idecal+k) = zr(iacocf-1+idecal_0+k)
                end do
                idecal = idecal+zi(iaconb-1+iNode2)
                cycle
            end if
        end if
        call pj2dap(iNode2, zr(iacoo2), zr(iacoo1), zi(iatr3), &
                    cobary, cellLink1, nbtrou, bt2ddi, bt2dvr, &
                    bt2dnb, bt2dlc, zi(iabtco), &
                    l_dmax, dmax, dala, loin, dmin)
        if (l_dmax .and. (nbtrou .eq. 0)) then
            zi(iaconb-1+iNode2) = 3
            zi(iacotr-1+iNode2) = 0
            cycle
        end if
        if (nbtrou .eq. 0) then
            nodeName2 = int_to_char8(iNode2)
            call utmess('F', 'PROJECTION4_56', sk=nodeName2)
        end if
        zi(iaconb-1+iNode2) = 3
        zi(iacotr-1+iNode2) = cellLink1
        do k = 1, 3
            zi(iaconu-1+idecal+k) = zi(iatr3+4*(cellLink1-1)+k)
            zr(iacocf-1+idecal+k) = cobary(k)
        end do
        idecal = idecal+zi(iaconb-1+iNode2)
    end do

! - For alarm, see ticket 16186

! - Transform corrMeshTemp in corrMesh (real cells)
    if (present(listInterc_)) then
        ASSERT(present(nbInterc_))
        listInterc = listInterc_
        nbInterc = nbInterc_
        call pj2dtr(corrMeshTemp, corrMesh, &
                    cellListType, cellListCode, &
                    zr(iacoo1), zr(iacoo2), &
                    spacedim, dala, &
                    listInterc, nbInterc)
    else
        listInterc = ' '
        nbInterc = 0
        call pj2dtr(corrMeshTemp, corrMesh, &
                    cellListType, cellListCode, &
                    zr(iacoo1), zr(iacoo2), &
                    spacedim, dala)
    end if

! - Debug
    if (dbg) then
        call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, '&&PJ2DCO', 1, ' ')
        call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, corrMesh, 1, ' ')
    end if

! - Clean
    call detrsd('CORRESP_2_MAILLA', corrMeshTemp)
    call jedetr(boite//'.BT2DDI')
    call jedetr(boite//'.BT2DVR')
    call jedetr(boite//'.BT2DNB')
    call jedetr(boite//'.BT2DLC')
    call jedetr(boite//'.BT2DCO')
    call jedetr('&&PJXXCO.LIMA1')
    call jedetr('&&PJXXCO.LIMA2')
    call jedetr('&&PJXXCO.LINO1')
    call jedetr('&&PJXXCO.LINO2')
    call jedetr('&&PJXXCO.TRIA3')
!
    call jedema()
end subroutine
