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
! ==================================================================================================
!
! Module for elementary computation of loads
!
! ==================================================================================================
!
module loadElemCompute_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: evalPresAtNodes, evalLineLoad
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! evalPresAtNodes
!
! Evaluation of pressure at nodes
!
! In  typeScal         : type of scalar for load (real, complex or function)
! In  nbNode           : number of nodes of element
! In  npg              : number of Gauss point
! In  jvPres           : adress to data from load
! In  jvTime           : adress to time parameters
! In  jvGeom           : adress to initial coordinates of nodes
! In  lAbsCurv         : flag for curvilinear coordinates
! In  absCurv          : value of curvilinear coordinates at nodes
! Out presNode         : value of pressure at nodes
! --------------------------------------------------------------------------------------------------
    subroutine evalPresAtNodes(typeScal, &
                               nbNode, &
                               jvPres, jvTime, jvGeom, &
                               lAbsCurv, absCurv, &
                               presNode)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: typeScal
        integer(kind=8), intent(in) :: nbNode
        integer(kind=8), intent(in) :: jvPres, jvTime, jvGeom
        aster_logical, intent(in) :: lAbsCurv
        real(kind=8), intent(in) :: absCurv(nbNode)
        real(kind=8), intent(out) :: presNode(nbNode)
! ----- Local
        integer(kind=8) :: iNode, ier
        integer(kind=8) :: nbPara
        integer(kind=8), parameter :: nbParaMax = 5
        character(len=8), parameter :: paraName(nbParaMax) = (/'X   ', 'Y   ', 'Z   ', &
                                                               'INST', 'ABSC'/)
        real(kind=8) :: paraVale(nbParaMax)
!   ------------------------------------------------------------------------------------------------
!
        presNode = 0.d0

! ----- Value of load at nodes
        if (typeScal .eq. 'R') then
            ASSERT(jvTime .eq. 0)
            do iNode = 1, nbNode
                presNode(iNode) = zr(jvPres-1+iNode)
            end do

        else if (typeScal .eq. 'F') then
            paraVale = 0.d0
            paraVale(4) = zr(jvTime)
            if (lAbsCurv) then
                nbPara = 5
            else
                nbPara = 4
            end if
            do iNode = 1, nbNode
                paraVale(1) = zr(jvGeom+3*(iNode-1))
                paraVale(2) = zr(jvGeom+3*(iNode-1)+1)
                paraVale(3) = zr(jvGeom+3*(iNode-1)+2)
                if (lAbsCurv) then
                    paraVale(5) = absCurv(iNode)
                end if
                call fointe('FM', zk8(jvPres), nbPara, paraName, paraVale, &
                            presNode(iNode), ier)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! evalLineLoad
!
! Evaluation of lineic force
!
! In  typeScal         : type of scalar for load (real, complex or function)
! In  lCplxRealPart    : if load if complex, flag for real part
! In  lGravity         : flag when lineic force is gravity
! In  iNode            : index of current node
! In  npg              : number of Gauss point
! In  jvLoad           : adress to data from load
! In  jvTime           : adress to time parameters
! In  jvGeom           : adress to initial coordinates of nodes
! Out lineLoad         : value of load for components DX, DY, DZ, DRX, DRY, DRZ
! In  rho              : density
! --------------------------------------------------------------------------------------------------
    subroutine evalLineLoad(typeScal, lCplxRealPart, lGravity, &
                            jvLoad, jvTime, jvGeom, &
                            rho, lineLoad, iNode_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: typeScal
        aster_logical, intent(in) :: lCplxRealPart, lGravity
        integer(kind=8), intent(in) :: jvTime, jvGeom, jvLoad
        real(kind=8), intent(in) :: rho
        real(kind=8), intent(out) :: lineLoad(6)
        integer(kind=8), optional, intent(in) :: iNode_
! ----- Local
        real(kind=8) :: gravAcce
        integer(kind=8), parameter :: nbParaMaxi = 4
        character(len=16), parameter :: paraName(nbParaMaxi) = (/'X   ', 'Y   ', 'Z   ', 'INST'/)
        real(kind=8) :: paraVale(nbParaMaxi)
        integer(kind=8) :: iret, nbPara, iDof
        character(len=8) :: funcName
!   ------------------------------------------------------------------------------------------------
!
        lineLoad = 0.d0
        if (lGravity) then
            ASSERT(lCplxRealPart)
            ASSERT(typeScal .eq. 'R')
            ASSERT(rho .gt. 0.d0)
            gravAcce = zr(jvLoad)
            lineLoad(1) = rho*gravAcce*zr(jvLoad+1)
            lineLoad(2) = rho*gravAcce*zr(jvLoad+2)
            lineLoad(3) = rho*gravAcce*zr(jvLoad+3)
        else
            if (typeScal .eq. 'R') then
                ASSERT(lCplxRealPart)
                lineLoad(1) = zr(jvLoad-1+1)
                lineLoad(2) = zr(jvLoad-1+2)
                lineLoad(3) = zr(jvLoad-1+3)
            else if (typeScal .eq. 'C') then
                lineLoad(1) = zr(jvLoad-1+1)
                lineLoad(2) = zr(jvLoad-1+2)
                lineLoad(3) = zr(jvLoad-1+3)
                if (lCplxRealPart) then
                    lineLoad(1) = dble(zc(jvLoad-1+1))
                    lineLoad(2) = dble(zc(jvLoad-1+2))
                    lineLoad(3) = dble(zc(jvLoad-1+3))
                else
                    lineLoad(1) = dimag(zc(jvLoad-1+1))
                    lineLoad(2) = dimag(zc(jvLoad-1+2))
                    lineLoad(3) = dimag(zc(jvLoad-1+3))
                end if
            elseif (typeScal .eq. 'F') then
                nbPara = 3
                if (jvTime .ne. 0) then
                    paraVale(4) = zr(jvTime)
                    nbPara = 4
                end if
                paraVale(1) = zr(jvGeom-1+3*(iNode_-1)+1)
                paraVale(2) = zr(jvGeom-1+3*(iNode_-1)+2)
                paraVale(3) = zr(jvGeom-1+3*(iNode_-1)+3)
                do iDof = 1, 6
                    funcName = zk8(jvLoad+iDof-1)
                    call fointe('FM', funcName, nbPara, paraName, paraVale, &
                                lineLoad(iDof), iret)
                end do
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module loadElemCompute_module
