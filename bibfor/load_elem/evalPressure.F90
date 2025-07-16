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
subroutine evalPressure(lFunc, lTime, time, &
                        nbNode, cellDime, ipg, &
                        jvShapFunc, jvGeom, jvPres, &
                        pres, cisa_, geomCurr_)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/fointe.h"
#include "asterfort/evalPressureSetFuncPara.h"
!
    aster_logical, intent(in) :: lFunc, lTime
    integer(kind=8), intent(in) :: cellDime, nbNode, ipg
    integer(kind=8), intent(in) :: jvGeom, jvShapFunc, jvPres
    real(kind=8), intent(in) :: time
    real(kind=8), intent(out) :: pres
    real(kind=8), optional, intent(out) :: cisa_
    real(kind=8), optional, intent(in) :: geomCurr_(*)
!
! --------------------------------------------------------------------------------------------------
!
! Load - Compute vector - PRES_REP
!
! Evaluation of pressure at Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  lFunc           : flag if load is function
! In  lTime           : flag if have time for function
! In  time            : value of current time
! In  nbNode          : total number of nodes
! In  cellDime        : dimension of cell (2 or 3)
! In  ipg             : current index of Gauss point
! In  jvShapFunc      : JEVEUX adress for shape functions
! In  jvGeom          : JEVEUX adress for geometry (coordinates of nodes)
! In  jvPres          : JEVEUX adress for pressure
! Out pres            : normal pressure
! Out cisa            : "tangent" pressure (shear, only in 2D !)
! In  geomCurr_       : cuurent geometry
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: paraNbMax = 7
    character(len=8) :: paraName(paraNbMax)
    real(kind=8) :: paraVale(paraNbMax)
    integer(kind=8) :: iNode, ldec, iret, paraNb
!
! --------------------------------------------------------------------------------------------------
!
    pres = 0.d0
    if (present(cisa_)) then
        cisa_ = 0.d0
    end if
!
    ldec = (ipg-1)*nbNode
!
    if (lFunc) then

! ----- Prepare parameters when pressure is function
        call evalPressureSetFuncPara(lTime, time, &
                                     nbNode, cellDime, ipg, &
                                     jvShapFunc, jvGeom, &
                                     paraNbMax, paraNb, paraName, paraVale, &
                                     geomCurr_)

! ----- Evaluate function
        call fointe('FM', zk8(jvPres-1+1), paraNb, paraName, paraVale, pres, iret)
        if (present(cisa_)) then
            call fointe('FM', zk8(jvPres-1+2), paraNb, paraName, paraVale, cisa_, iret)
        end if
    else
        if (present(cisa_)) then
            do iNode = 1, nbNode
                pres = pres+zr(jvPres+2*(iNode-1)-1+1)*zr(jvShapFunc+ldec-1+iNode)
                cisa_ = cisa_+zr(jvPres+2*(iNode-1)-1+2)*zr(jvShapFunc+ldec-1+iNode)
            end do
        else
            do iNode = 1, nbNode
                pres = pres+zr(jvPres+(iNode-1)-1+1)*zr(jvShapFunc+ldec-1+iNode)
            end do
        end if
    end if
!
end subroutine
