! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine pjNodeNode(meshZ, iOcc, coniJvZ,&
                      tran, cent, rotaMatr,&
                      nbNodeMast, nbNodeSlav,&
                      nodeMast, nodeSlav)
!
implicit none
!
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/padist.h"
#include "asterfort/parotr.h"
#include "asterfort/utmess.h"
!
character(len=*), intent(in) :: meshZ
integer, intent(in) :: iOcc
character(len=*), intent(in) :: coniJvZ
real(kind=8), intent(in) :: tran(3), cent(3), rotaMatr(3, 3)
integer, intent(in) :: nbNodeMast, nbNodeSlav
integer, pointer :: nodeMast(:), nodeSlav(:)
!
! --------------------------------------------------------------------------------------------------
!
! Compute projection of node on another node
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: lDebug = ASTER_FALSE
    character(len=8) :: mesh, nameNodeMast, nameNodeSlav
    character(len=24) :: coniJv
    integer :: nbNode, iNodeMast, iNodeSlav, iNode1, iNode2, iNode
    integer :: numeNodeMast, numeNodeSlav, linkIndx, numeNodeLink, nbError
    integer :: jvGeom
    integer, pointer :: coni(:) => null()
    real(kind=8) :: coorMast(3), coorSlav(3)
    real(kind=8) :: distMini, dist
    aster_logical :: lNodeLinked
    integer, pointer :: nodeInve(:) => null()
    integer, pointer :: nodeOut1(:) => null(), nodeOut2(:) => null()
    integer, pointer :: nodeOut3(:) => null(), nodeOut4(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    mesh = meshZ
    coniJv = coniJvZ

! - Access to mesh
    call jeveuo(mesh//'.COORDO    .VALE', 'L', jvGeom)

! - Number of nodes
    ASSERT(nbNodeMast .eq. nbNodeSlav)
    nbNode = nbNodeMast

! - Working objects
    AS_ALLOCATE(vi = nodeInve, size = nbNode)
    AS_ALLOCATE(vi = nodeOut1, size = nbNode)
    AS_ALLOCATE(vi = nodeOut2, size = nbNode)
    AS_ALLOCATE(vi = nodeOut3, size = nbNode)
    AS_ALLOCATE(vi = nodeOut4, size = nbNode)

! - First link: from master to slave
    nodeInve = 0
    nbError = 0
    do iNodeMast = 1, nbNode
! ----- Current master node
        numeNodeMast = nodeMast(iNodeMast)

! ----- Apply transformation on master node
        call parotr(mesh, jvGeom, numeNodeMast, 0, cent,&
                    rotaMatr, tran, coorMast)

! ----- Find slave node coupled to master node
        distMini = r8gaem()
        linkIndx = 0
        do iNodeSlav = 1, nbNode
! --------- Current slave node
            numeNodeSlav = nodeSlav(iNodeSlav)
            coorSlav(1) = zr(jvGeom-1+3*(numeNodeSlav-1)+1)
            coorSlav(2) = zr(jvGeom-1+3*(numeNodeSlav-1)+2)
            coorSlav(3) = zr(jvGeom-1+3*(numeNodeSlav-1)+3)

! --------- Distance and test
            dist = padist( 3, coorMast, coorSlav )
            if (dist .lt. distMini) then
                distMini = dist
                linkIndx = iNodeSlav
                numeNodeLink = numeNodeSlav
            endif
        end do

! ----- No node found
        if (linkIndx .eq. 0) then
            call jenuno(jexnum(mesh//'.NOMNOE', numeNodeMast), nameNodeMast)
            call utmess('F', 'CHARGES7_3', sk=nameNodeMast)
        endif

! ----- Current node
        if (nodeInve(linkIndx) .eq. 0) then
            nodeOut1(iNodeMast) = numeNodeMast
            nodeOut2(iNodeMast) = numeNodeLink
            nodeInve(linkIndx) = numeNodeMast
        else
            nbError = nbError + 1
            call utmess('E', 'CHARGES7_77')
        endif
    end do

! - Debug
    if (lDebug) then
        WRITE(6,*) 'First link: from master to slave'
        do iNode = 1, nbNode
            WRITE(6,*) "Node ",nodeOut1(iNode),' with ',nodeOut2(iNode)
        end do
    endif

! - Final error
    if (nbError .ne. 0) then
        call utmess('F', 'CHARGES7_7')
    endif

! - Second link: from slave to master
    nodeInve = 0
    nbError = 0
    do iNodeSlav = 1, nbNode
! ----- Current slave node
        numeNodeSlav = nodeSlav(iNodeSlav)
        coorSlav(1) = zr(jvGeom-1+3*(numeNodeSlav-1)+1)
        coorSlav(2) = zr(jvGeom-1+3*(numeNodeSlav-1)+2)
        coorSlav(3) = zr(jvGeom-1+3*(numeNodeSlav-1)+3)

! ----- Find master node coupled to slave node
        distMini = r8gaem()
        linkIndx = 0
        do iNodeMast = 1, nbNode
            numeNodeMast = nodeMast(iNodeMast)

! --------- Apply transformation on master node
            call parotr(mesh, jvGeom, numeNodeMast, 0, cent,&
                        rotaMatr, tran, coorMast)

! --------- Distance and test
            dist = padist( 3, coorMast, coorSlav )
            if (dist .lt. distMini) then
                distMini = dist
                linkIndx = iNodeMast
                numeNodeLink = numeNodeMast
            endif
        end do

! ----- No node found
        if (linkIndx .eq. 0) then
            call jenuno(jexnum(mesh//'.NOMNOE', numeNodeSlav), nameNodeSlav)
            call utmess('F', 'CHARGES7_3', sk=nameNodeSlav)
        endif

! ----- Current node
        if (nodeInve(linkIndx) .eq. 0) then
            nodeOut3(iNodeSlav) = numeNodeLink
            nodeOut4(iNodeSlav) = numeNodeSlav
            nodeInve(linkIndx) = numeNodeLink
        else
            nbError = nbError + 1
            call utmess('E', 'CHARGES7_77')
        endif
    end do

! - Final error
    if (nbError .ne. 0) then
        call utmess('F', 'CHARGES7_7')
    endif

! - Debug
    if (lDebug) then
        WRITE(6,*) 'Second link: from slave to master'
        do iNode = 1, nbNode
            WRITE(6,*) "Node ",nodeOut3(iNode),' with ',nodeOut4(iNode)
        end do
    endif

! - Check consistency of the two lists
    do iNode1 = 1, nbNode
        lNodeLinked = ASTER_FALSE
        do iNode2 = 1, nbNode
            if (nodeOut1(iNode1) .eq. nodeOut3(iNode2)) then
                lNodeLinked = ASTER_TRUE
                if (nodeOut2(iNode1) .ne. nodeOut4(iNode2)) then
                    nbError = nbError + 1
                    call utmess('E', 'CHARGES7_88')
                endif
            endif
        end do
        if (.not. lNodeLinked) then
            nbError = nbError + 1
            call utmess('E', 'CHARGES7_88')
        endif
    end do

! - Final error
    if (nbError .ne. 0) then
        call utmess('F', 'CHARGES7_7')
    endif

! - Write final pairing
    call jecroc(jexnum(coniJv, iOcc))
    call jeecra(jexnum(coniJv, iOcc), 'LONMAX', 2*nbNode+1)
    call jeveuo(jexnum(coniJv, iOcc), 'E', vi = coni)
    coni(1) = nbNode
    do iNode = 1, nbNode
        coni(1+2*(iNode-1)+1) = nodeOut1(iNode)
        coni(1+2*(iNode-1)+2) = nodeOut2(iNode)
    end do

! - Clean
    AS_DEALLOCATE(vi = nodeInve)
    AS_DEALLOCATE(vi = nodeOut1)
    AS_DEALLOCATE(vi = nodeOut2)
    AS_DEALLOCATE(vi = nodeOut3)
    AS_DEALLOCATE(vi = nodeOut4)
end subroutine
