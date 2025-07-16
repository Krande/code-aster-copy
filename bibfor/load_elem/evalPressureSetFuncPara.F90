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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine evalPressureSetFuncPara(lTime, time, &
                                   nbNode, cellDime, ipg, &
                                   jvShapFunc, jvGeom, &
                                   paraNbMax, paraNb, paraName, paraVale, &
                                   geomCurr_)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
!
    aster_logical, intent(in) :: lTime
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: nbNode, cellDime, ipg
    integer(kind=8), intent(in) :: jvShapFunc, jvGeom
    integer(kind=8), intent(in) :: paraNbMax
    integer(kind=8), intent(out) :: paraNb
    character(len=8) :: paraName(paraNbMax)
    real(kind=8) :: paraVale(paraNbMax)
    real(kind=8), optional, intent(in) :: geomCurr_(*)
!
! --------------------------------------------------------------------------------------------------
!
! Load - Compute vector - PRES_REP
!
! Prepare parameters when pressure is function
!
! --------------------------------------------------------------------------------------------------
!
! In  lTime           : flag if have time for function
! In  time            : value of current time
! In  nbNode          : total number of nodes
! In  cellDime        : dimension of cell (2 or 3)
! In  ipg             : current index of Gauss point
! In  jvShapFunc      : JEVEUX adress for shape functions
! In  jvGeom          : JEVEUX adress for geometry (coordinates of nodes)
! In  geom_reac       : updated geometry
! Out paraNb          : number of parameters to evaluate function
! Out paraName        : name of parameters to evaluate function
! Out paraVale        : value of parameters to evaluate function
! In  geomCurr_       : cuurent geometry
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: x, y, z, xf, yf, zf
    integer(kind=8) :: iNode, ldec
!
! --------------------------------------------------------------------------------------------------
!
    ldec = (ipg-1)*nbNode

! - Initial and current coordinates of current Gauss point
    x = 0.d0
    y = 0.d0
    z = 0.d0
    xf = 0.d0
    yf = 0.d0
    zf = 0.d0
    do iNode = 1, nbNode
        x = x+zr(jvGeom+(cellDime+1)*(iNode-1)-1+1)*zr(jvShapFunc+ldec-1+iNode)
        y = y+zr(jvGeom+(cellDime+1)*(iNode-1)-1+2)*zr(jvShapFunc+ldec-1+iNode)
        if (present(geomCurr_)) then
            xf = xf+geomCurr_((cellDime+1)*(iNode-1)+1)*zr(jvShapFunc+ldec-1+iNode)
            yf = yf+geomCurr_((cellDime+1)*(iNode-1)+2)*zr(jvShapFunc+ldec-1+iNode)
        end if
        if (cellDime .eq. 2) then
            z = z+zr(jvGeom+(cellDime+1)*(iNode-1)-1+3)*zr(jvShapFunc+ldec-1+iNode)
            if (present(geomCurr_)) then
                zf = zf+geomCurr_((cellDime+1)*(iNode-1)+3)*zr(jvShapFunc+ldec-1+iNode)
            end if
        end if
    end do

! - List of parameters
    paraNb = 4
    paraVale(1) = x
    paraName(1) = 'X'
    paraVale(2) = xf
    paraName(2) = 'XF'
    paraVale(3) = y
    paraName(3) = 'Y'
    paraVale(4) = yf
    paraName(4) = 'YF'
    if (cellDime .eq. 2) then
        paraNb = paraNb+1
        paraVale(paraNb) = z
        paraName(paraNb) = 'Z'
        paraNb = paraNb+1
        paraVale(paraNb) = zf
        paraName(paraNb) = 'ZF'
    end if
    if (lTime) then
        paraNb = paraNb+1
        paraVale(paraNb) = time
        paraName(paraNb) = 'INST'
    end if
!
end subroutine
