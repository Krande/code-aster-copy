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
subroutine te0384(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/evalFaceSpeedDire.h"
#include "asterfort/evalFaceSpeedVale.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vff2dn.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D_FLUI_STRU, AXIS_FLUI_STRU
!
! Options: CHAR_MECA_VFAC
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lFunc, lTime
    integer(kind=8) :: jvGeom, jvMate, jvLoad, jvTime, jvVect
    aster_logical :: lReal
    real(kind=8) :: x, y, speedDire
    real(kind=8) :: nx, ny
    real(kind=8) :: rho, poids
    real(kind=8) :: time, speedVale
    integer(kind=8) :: jvWeight, jvShape, jvDShape
    integer(kind=8) :: nbNode, npg, cellDime, ndofbynode
    integer(kind=8) :: ldec
    integer(kind=8) :: i, ii, ipg
    aster_logical :: l_axis
    real(kind=8) :: r
    integer(kind=8) :: j_mater, iret
    character(len=16) :: FEForm
!
! --------------------------------------------------------------------------------------------------
!
    lFunc = (option .eq. 'CHAR_MECA_VFAC_F')
    lReal = .not. lFunc

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMate)
    if (lFunc) then
        call jevech('PVITEFF', 'L', jvLoad)
    else
        call jevech('PVITEFR', 'L', jvLoad)
    end if

! - Get time if present
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    lTime = ASTER_FALSE
    time = 0.d0
    if (jvTime .ne. 0) then
        lTime = ASTER_TRUE
        time = zr(jvTime)
    end if

! - Get element parameters
    l_axis = (lteatt('AXIS', 'OUI'))
    call teattr('S', 'FORMULATION', FEForm, iret)
    call elrefe_info(fami='RIGI', &
                     nno=nbNode, npg=npg, ndim=cellDime, &
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShape)
    ASSERT(nbNode .le. 3)
    if (FEForm .eq. 'U_P_PHI') then
        ndofbynode = 3
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        ndofbynode = 3
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if

! - Get material properties for fluid
    j_mater = zi(jvMate)
    call getFluidPara(j_mater, rho_=rho)

! - Output field
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, ndofbynode
        zr(jvVect+i-1) = 0.d0
    end do

! - Loop on Gauss points
    do ipg = 1, npg
        ldec = (ipg-1)*nbNode

! ----- Compute normal
        nx = 0.d0
        ny = 0.d0
        call vff2dn(cellDime, nbNode, ipg, jvWeight, jvDShape, &
                    zr(jvGeom), nx, ny, poids)
        if (l_axis) then
            r = 0.d0
            do i = 1, nbNode
                r = r+zr(jvGeom+2*(i-1))*zr(jvShape+ldec+i-1)
            end do
            poids = poids*r
        end if

! ----- Get value of speed
        call evalFaceSpeedVale(lFunc, lTime, time, &
                               nbNode, cellDime, ipg, &
                               jvShape, jvGeom, jvLoad, &
                               speedVale, x, y)

! ----- Get direction of speed
        call evalFaceSpeedDire(FEForm, cellDime, jvLoad, speedDire, &
                               nx, ny, &
                               lFunc_=lFunc, lReal_=lReal, &
                               lTime_=lTime, time_=time, &
                               x_=x, y_=y)

! ----- Normal convention is fluid -> structure
        speedDire = -speedDire

! ----- Compute vector
        do i = 1, nbNode
            ii = ndofbynode*i
            zr(jvVect+ii-1) = zr(jvVect+ii-1)-speedDire* &
                              poids* &
                              zr(jvShape+ldec+i-1)*speedVale*rho
        end do

    end do
!
end subroutine
