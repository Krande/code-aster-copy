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
subroutine cmmoma(meshOutZ, nbCellModi, modiCellNume, modiCellType, nbNodeIn)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!
    character(len=*), intent(in) :: meshOutZ
    integer(kind=8), intent(in) :: nbNodeIn, nbCellModi
    integer(kind=8), pointer :: modiCellNume(:), modiCellType(:)
!
!     OPERATEUR CREA_MAILLAGE   MOT CLE FACTEUR "MODI_MAILLE"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: cellNume, cellType, nodeAddIndx, nodeNume
    integer(kind=8) :: iDime, iCell
    real(kind=8) :: w
    character(len=8) :: meshOut
    real(kind=8) :: coo1(3), coo2(3), coo3(3), theta, epsi, t13(3), t32(3)
    real(kind=8) :: normen, norme1, norme2, n(3), om(3), oc(3), c2, c6, t2, t6
    real(kind=8) :: t12(3)
    real(kind=8) :: n3m(3), mc(3), mp(3), mr(3), x3(3), x4(3), costet, dn1n2
    integer(kind=8) :: icoude, iNode, nodeNume1, nodeNume2, nodeNume3, nodeNume4, cellModiType
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
!
    meshOut = meshOutZ
!
! - Access to mesh
!
    call jeveuo(meshOut//'.TYPMAIL', 'E', vi=typmail)
    call jeveuo(meshOut//'.COORDO    .VALE', 'E', vr=vale)
!
    do iCell = 1, nbCellModi
! ----- Cell to modify
        cellNume = modiCellNume(iCell)
        cellType = modiCellType(iCell)
!
! ----- New type of cell
        if (cellType .eq. MT_TRIA6) then
            cellModiType = MT_TRIA7
            nodeAddIndx = 7
        else if (cellType .eq. MT_QUAD8) then
            cellModiType = MT_QUAD9
            nodeAddIndx = 9
        else if (cellType .eq. MT_SEG3) then
            cellModiType = MT_SEG4
            nodeAddIndx = 4
        else
            ASSERT(ASTER_FALSE)
        end if
        typmail(cellNume) = cellModiType
!
        call jeveuo(jexnum(meshOut//'.CONNEX', cellNume), 'E', vi=connex)
!
! ----- Add node
        nodeNume = nbNodeIn+iCell
        connex(nodeAddIndx) = nodeNume
!
        if (cellType .eq. MT_SEG3) then
            nodeNume4 = nodeNume
! --------- Coordinates of nodes
            do iNode = 1, 3
                nodeNume1 = connex(1)
                coo1(iNode) = vale(3*(nodeNume1-1)+iNode)
                nodeNume2 = connex(2)
                coo2(iNode) = vale(3*(nodeNume2-1)+iNode)
                nodeNume3 = connex(3)
                coo3(iNode) = vale(3*(nodeNume3-1)+iNode)
            end do
!
            t13 = coo3-coo1
            t32 = coo2-coo3
            t12 = coo2-coo1
            call normev(t13, norme1)
            call normev(t32, norme2)
            call normev(t12, dn1n2)
            call provec(t32, t13, n)
            call normev(n, normen)
            epsi = 1.d-4*norme1
!
!           VERIF QUE LE 3EME NOEUD EST BIEN AU MILIEU
!
            if (abs(norme2-norme1) .gt. epsi) then
                call utmess('F', 'MESH2_23')
            end if
!
            if (normen .le. epsi) then
                icoude = 0
                theta = 0.d0
            else
                icoude = 1
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                costet = ddot(b_n, t13, b_incx, t32, b_incy)
                theta = 2.d0*atan2(normen, costet)
            end if
!
            if (icoude .eq. 0) then
                do iNode = 1, 3
                    x3(iNode) = coo1(iNode)+t12(iNode)*dn1n2/3.d0
                    x4(iNode) = coo1(iNode)+2.d0*t12(iNode)*dn1n2/3.d0
                end do
            else
                c2 = cos(theta/2.d0)
                c6 = cos(theta/6.d0)
                t2 = tan(theta/2.d0)
                t6 = tan(theta/6.d0)
                do iNode = 1, 3
                    om(iNode) = (coo1(iNode)+coo2(iNode))*0.5d0
                    n3m(iNode) = om(iNode)-coo3(iNode)
                    mc(iNode) = n3m(iNode)*c2/(1.d0-c2)
                    oc(iNode) = om(iNode)+mc(iNode)
                    mp(iNode) = (coo1(iNode)-om(iNode))*t6/t2
                    mr(iNode) = (coo2(iNode)-om(iNode))*t6/t2
                    x3(iNode) = oc(iNode)+(mp(iNode)-mc(iNode))*c6/c2
                    x4(iNode) = oc(iNode)+(mr(iNode)-mc(iNode))*c6/c2
                end do
            end if
            vale(3*(nodeNume3-1)+1) = x3(1)
            vale(3*(nodeNume3-1)+2) = x3(2)
            vale(3*(nodeNume3-1)+3) = x3(3)
            vale(3*(nodeNume4-1)+1) = x4(1)
            vale(3*(nodeNume4-1)+2) = x4(2)
            vale(3*(nodeNume4-1)+3) = x4(3)
!
        else if (cellType .eq. MT_TRIA6) then
            do iDime = 1, 3
                w = 0.d0
                w = w+vale(3*(connex(1)-1)+iDime)*(-1.d0/9.d0)
                w = w+vale(3*(connex(2)-1)+iDime)*(-1.d0/9.d0)
                w = w+vale(3*(connex(3)-1)+iDime)*(-1.d0/9.d0)
                w = w+vale(3*(connex(4)-1)+iDime)*(4.d0/9.d0)
                w = w+vale(3*(connex(5)-1)+iDime)*(4.d0/9.d0)
                w = w+vale(3*(connex(6)-1)+iDime)*(4.d0/9.d0)
                vale(3*(nodeNume-1)+iDime) = w
            end do
        else if (cellType .eq. MT_QUAD8) then
            do iDime = 1, 3
                w = 0.d0
                w = w+vale(3*(connex(1)-1)+iDime)*(-1.d0/4.d0)
                w = w+vale(3*(connex(2)-1)+iDime)*(-1.d0/4.d0)
                w = w+vale(3*(connex(3)-1)+iDime)*(-1.d0/4.d0)
                w = w+vale(3*(connex(4)-1)+iDime)*(-1.d0/4.d0)
                w = w+vale(3*(connex(5)-1)+iDime)*(1.d0/2.d0)
                w = w+vale(3*(connex(6)-1)+iDime)*(1.d0/2.d0)
                w = w+vale(3*(connex(7)-1)+iDime)*(1.d0/2.d0)
                w = w+vale(3*(connex(8)-1)+iDime)*(1.d0/2.d0)
                vale(3*(nodeNume-1)+iDime) = w
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
end subroutine
