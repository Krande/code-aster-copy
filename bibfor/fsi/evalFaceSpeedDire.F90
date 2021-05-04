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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine evalFaceSpeedDire(fsi_form, cellDime, jvLoad , speedDire,&
                             ipg     , nx     , ny       ,&
                             lFunc_  , lReal_ , lCplx_   ,&
                             lTime_  , time_  ,&
                             x_      , y_     ,&
                             z_      , nz_)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/utmess.h"
!
character(len=16), intent(in) :: fsi_form
integer, intent(in) :: cellDime, jvLoad
real(kind=8), intent(out) :: speedDire
integer, intent(in) :: ipg
real(kind=8), intent(in) :: nx, ny
aster_logical, optional, intent(in) :: lFunc_, lReal_, lCplx_, lTime_
real(kind=8), optional, intent(in) :: time_, x_, y_
real(kind=8), optional, intent(in) :: z_, nz_
!
! --------------------------------------------------------------------------------------------------
!
! Utilities for fluid
!
! Evaluation of direction of speed for VITE_FACE
!
! --------------------------------------------------------------------------------------------------
!
! In  fsi_form         : type of FORMULATION
! In  cellDime         : dimension of cell (2 or 3)
! In  jvLoad           : JEVEUX adress for field with parameters for load
! Out speedDire        : direction of speed (dot product with normal)
! In  lFunc            : flag if VITE_FACE is function
! In  lReal            : flag if VITE_FACE is real
! In  lCplx            : flag if VITE_FACE is complex
! In  x, y, z          : coordinates of current Gauss point
! In  lTime            : flag if have time for function
! In  time             : value of current time
! In  nx, ny, nz       : normal to face
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbCmpInMap = 5
    aster_logical :: lSpeedNorm
    aster_logical :: lFunc, lReal, lCplx, lTime
    character(len=8) :: funcName
    real(kind=8) :: time, x, y, z, nz
    real(kind=8) :: nvx, nvy, nvz, normVect
    real(kind=8) :: nxNorm, nyNorm, nzNorm
    integer :: nbPara, iret
    character(len=8), parameter :: paraName(4) = (/'X   ', 'Y   ',  'Z   ', 'INST'/)
    real(kind=8) :: paraVale(4)
!
! --------------------------------------------------------------------------------------------------
!
    lFunc = ASTER_FALSE
    lReal = ASTER_FALSE
    lCplx = ASTER_FALSE
    lTime = ASTER_FALSE
    if (present(lFunc_)) lFunc = lFunc_
    if (present(lReal_)) lReal = lReal_
    if (present(lCplx_)) lCplx = lCplx_
    if (present(lTime_)) lTime = lTime_
    x = 0.d0
    y = 0.d0
    z = 0.d0
    if (present(x_)) x = x_
    if (present(y_)) y = y_
    if (present(z_)) z = z_
    nz = 0.d0
    if (present(nz_)) nz = nz_
    time = 0.d0
    if (present(time_)) time = time_
    speedDire = 1.d0

! - Norm vector
    normVect = sqrt(nx*nx + ny*ny + nz*nz)
    if (normVect .le. r8prem()) then
        call utmess('F', 'CHARGES6_6')
    endif
    nxNorm = nx / normVect
    nyNorm = ny / normVect
    nzNorm = nz / normVect

! - Flag if normal speed or not
    lSpeedNorm = ASTER_TRUE
    if (lFunc) then
        lSpeedNorm = zk8(jvLoad-1+nbCmpInMap*(ipg-1)+2) .eq. '&FOZERO'

    elseif (lReal) then
        lSpeedNorm = zr(jvLoad-1+nbCmpInMap*(ipg-1)+2) .le. 0.d0

    elseif (lCplx) then
        lSpeedNorm = dble(zc(jvLoad-1+nbCmpInMap*(ipg-1)+2)) .le. 0.d0

    else
        ASSERT(ASTER_FALSE)

    endif

    if (lSpeedNorm) then
        speedDire = 1.d0

    else
        if (fsi_form .ne. 'FSI_UPSI') then
            call utmess('F', 'CHARGES6_7', sk = fsi_form)
        endif
        nvx = 0.d0
        nvy = 0.d0
        nvz = 0.d0
        if (lFunc) then
            nbPara = 2
            paraVale(1) = x
            paraVale(2) = y
            if (cellDime .eq. 2) then
                nbPara = 3
                paraVale(3) = z
            endif
            if (lTime) then
                nbPara = nbPara + 1
                paraVale(nbPara) = time
            endif
            funcName = zk8(jvLoad-1+nbCmpInMap*(ipg-1)+3)
            call fointe('FM', funcName, nbPara, paraname, paraVale, nvx, iret)
            funcName = zk8(jvLoad-1+nbCmpInMap*(ipg-1)+4)
            call fointe('FM', funcName, nbPara, paraname, paraVale, nvy, iret)
            if (cellDime .eq. 2) then
                funcName = zk8(jvLoad-1+nbCmpInMap*(ipg-1)+5)
                call fointe('FM', funcName, nbPara, paraname, paraVale, nvz, iret)
            endif

        elseif (lReal) then
            nvx = zr(jvLoad-1+nbCmpInMap*(ipg-1)+3)
            nvy = zr(jvLoad-1+nbCmpInMap*(ipg-1)+4)
            if (cellDime .eq. 2) then
                nvz = zr(jvLoad-1+nbCmpInMap*(ipg-1)+5)
            endif

        elseif (lCplx) then
            nvx = dble(zc(jvLoad-1+nbCmpInMap*(ipg-1)+3))
            nvy = dble(zc(jvLoad-1+nbCmpInMap*(ipg-1)+4))
            if (cellDime .eq. 2) then
                nvz = dble(zc(jvLoad-1+nbCmpInMap*(ipg-1)+5))
            endif

        endif

        normVect = sqrt(nvx*nvx + nvy*nvy + nvz*nvz)
        if (normVect .le. r8prem()) then
            speedDire = 1.d0
        else
            nvx = nvx / normVect
            nvy = nvy / normVect
            nvz = nvz / normVect
            if (cellDime .eq. 2) then
                speedDire = 1.d0*(nxNorm*nvx + nyNorm*nvy + nzNorm*nvz)
            else
                speedDire = 1.d0*(nxNorm*nvx + nyNorm*nvy )
            endif
        endif

    endif
!
end subroutine
