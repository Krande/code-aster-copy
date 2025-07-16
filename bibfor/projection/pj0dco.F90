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
subroutine pj0dco(typeSelect, &
                  entity1, entity2, &
                  nbCellSelect1, listCellSelect1, &
                  nbNodeSelect2, listNodeSelect2, &
                  geom1, geom2, corrMesh, &
                  dmax0d)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
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
    real(kind=8), intent(in) :: dmax0d

!
! --------------------------------------------------------------------------------------------------
!
! Create corresp_2_mailla datastructure
!
! Projection of nodes from entity1 to entity2 when entity1 is 0D
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
!    aster_logical, parameter :: dbg = ASTER_FALSE
    character(len=8) :: mesh1, mesh2, nodeName2, cdim1
    integer(kind=8) :: nbCellType
    integer(kind=8) :: cellListType(MT_NTYMAX)
    character(len=8) :: cellListCode(MT_NTYMAX)
    integer(kind=8) :: ifm, niv, nbNode2, nbCell1, k
    integer(kind=8) :: nbPoi, iase2, ndim, INode1, iPoi
    integer(kind=8) :: iCell1, iNode2
    integer(kind=8) :: iacoo1, iacoo2
    integer(kind=8) :: jxxk1, iaconu, iacocf, iacom1
    integer(kind=8) :: ilcnx1
    integer(kind=8) :: iaconb, cellTypeNume, idecal, cellLink1, nodeLink1
    real(kind=8) :: d2, dp, coor2(3), coor1(3), v12(3)
    integer(kind=8), pointer :: listCell1(:) => null(), listNode2(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

    ASSERT(typeSelect .eq. 'PARTIE')

! - Prepare list of entities

    call pjxxut('0D', typeSelect, &
                entity1, entity2, &
                nbCellSelect1, listCellSelect1, &
                nbNodeSelect2, listNodeSelect2, &
                mesh1, mesh2, &
                nbCellType, cellListType, cellListCode)

    call jeveuo('&&PJXXCO.LIMA1', 'L', vi=listCell1)
    call jeveuo('&&PJXXCO.LINO2', 'L', vi=listNode2)

! - Space dimension
    call dismoi('Z_QUASI_ZERO', mesh1, 'MAILLAGE', repk=cdim1)
    if (cdim1 .eq. 'OUI') then
        ndim = 2
    else
        ndim = 3
    end if

! - Access to meshes
    call dismoi('NB_NO_MAILLA', mesh2, 'MAILLAGE', repi=nbNode2)
    call dismoi('NB_MA_MAILLA', mesh1, 'MAILLAGE', repi=nbCell1)

! - Number of POI1
    call jeveuo(mesh1//'.TYPMAIL', 'L', vi=typmail)
    nbPoi = 0
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            cellTypeNume = typmail(iCell1)
            if (cellTypeNume .eq. cellListType(1)) then
                nbPoi = nbPoi+1
            else
                call utmess('F', 'PROJECTION4_1')
            end if
        end if
    end do
    if (nbPoi .eq. 0) then
        call utmess('F', 'PROJECTION4_55')
    end if

! - Create object of couples iCell1/Node1
!      on cree l'objet v='&&pjxxco.poi1' : ojb s v i
!         long(v)=1+2*nbPoi
!         v(1) : nbPoi(=nombre de poi1)
!         v(1+2(i-1)+1) : numero du noeud support
!         v(1+2(i-1)+2) : numero de la maille
!
    call wkvect('&&PJXXCO.POI1', 'V V I', 1+2*nbPoi, iase2)
    zi(iase2-1+1) = nbPoi
    call jeveuo(mesh1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    nbPoi = 0
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            cellTypeNume = typmail(iCell1)
            if (cellTypeNume .eq. cellListType(1)) then
                nbPoi = nbPoi+1
                zi(iase2+(nbPoi-1)*2+2) = iCell1
                zi(iase2+(nbPoi-1)*2+1) = connex(1+zi(ilcnx1-1+iCell1)-1)
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

! - Create temporary datastructure
    call wkvect(corrMesh//'.PJXX_K1', 'V V K24', 5, jxxk1)
    zk24(jxxk1-1+1) = mesh1
    zk24(jxxk1-1+2) = mesh2
    zk24(jxxk1-1+3) = 'COLLOCATION'
    call wkvect(corrMesh//'.PJEF_NB', 'V V I', nbNode2, iaconb)
    call wkvect(corrMesh//'.PJEF_NU', 'V V I', nbNode2, iaconu)
    call wkvect(corrMesh//'.PJEF_CF', 'V V R', nbNode2, iacocf)
    call wkvect(corrMesh//'.PJEF_M1', 'V V I', nbNode2, iacom1)

! - Look for cell cellLink1 in entity1 for each node from entity 2
    idecal = 0
    do iNode2 = 1, nbNode2
        if (listNode2(iNode2) .eq. 0) cycle
        coor2(1:ndim) = zr(iacoo2+ndim*(iNode2-1)-1+1:iacoo2+ndim*(iNode2-1)-1+ndim)
        d2 = r8maem()
        do iPoi = 1, nbPoi
            iNode1 = zi(iase2+(iPoi-1)*2+1)
            coor1(1:ndim) = zr(iacoo1+ndim*(iNode1-1)-1+1:iacoo1+ndim*(iNode1-1)-1+ndim)
            v12(1:ndim) = coor2(1:ndim)-coor1(1:ndim)
            dp = v12(1)**2+v12(2)**2
            if (ndim .eq. 3) dp = dp+v12(3)**2
            dp = sqrt(dp)
            if (dp .lt. d2) then
                d2 = dp
                cellLink1 = zi(iase2+(iPoi-1)*2+2)
                nodeLink1 = iNode1
            end if
        end do

        if (d2 .gt. dmax0d) then
            nodeName2 = int_to_char8(iNode2)
            call utmess('F', 'PROJECTION4_2', sk=nodeName2, nr=2, valr=[d2, dmax0d])
        end if

        zi(iaconb-1+iNode2) = 1
        zi(iacom1-1+iNode2) = cellLink1
        do k = 1, 2
            zi(iaconu-1+idecal+k) = nodeLink1
            zr(iacocf-1+idecal+k) = 1.d0
        end do
        idecal = idecal+1
    end do

! - Debug
!    if (dbg) then
!        call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, '&&PJ0DCO', 1, ' ')
!        call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, corrMesh, 1, ' ')
!    endif

    call jedetr('&&PJXXCO.LIMA1')
    call jedetr('&&PJXXCO.LIMA2')
    call jedetr('&&PJXXCO.LINO1')
    call jedetr('&&PJXXCO.LINO2')
    call jedetr('&&PJXXCO.POI1')
!
    call jedema()
end subroutine
