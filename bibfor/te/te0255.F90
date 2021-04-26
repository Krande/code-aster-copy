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
!
subroutine te0255(option, nomte)
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
! Elements: 2D_FLUIDE, AXIS_FLUIDE (boundary)
!
! Options: CHAR_MECA_VFAC
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lFunc, lTime
    integer :: jvGeom, jvMate, jvLoad, jvTime, jvVect
    aster_logical :: lReal
    real(kind=8) :: x, y, speedDire
    real(kind=8) :: nx, ny
    real(kind=8) :: rho, poids
    real(kind=8) :: time, speedVale
    integer :: jvWeight, jvShape, jvDShape
    integer :: nbNode, npg, cellDime, ndofbynode
    integer :: ldec
    integer :: i, ii, ipg
    aster_logical :: l_axis
    real(kind=8) :: r
    integer :: j_mater, iret
    character(len=16) :: fsi_form
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
    endif

! - Get time if present
    call tecach('NNO', 'PTEMPSR', 'L', iret, iad=jvTime)
    lTime = ASTER_FALSE
    time   = 0.d0
    if (jvTime .ne. 0) then
        lTime = ASTER_TRUE
        time   = zr(jvTime)
    endif

! - Get element parameters
    l_axis = (lteatt('AXIS','OUI'))
    call teattr('S', 'FORMULATION', fsi_form, iret)
    call elrefe_info(fami='RIGI',&
                     nno=nbNode, npg=npg, ndim=cellDime,&
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShape)
    ASSERT(nbNode .le. 3)
    if (fsi_form .eq. 'FSI_UPPHI') then
        ndofbynode = 2
    elseif (fsi_form .eq. 'FSI_UP') then
        ndofbynode = 1
    else
        call utmess('F', 'FLUID1_2', sk = fsi_form)
    endif

! - Get material properties for fluid
    j_mater = zi(jvMate)
    call getFluidPara(j_mater, rho)

! - Output field
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, ndofbynode*nbNode
        zr(jvVect+i-1) = 0.d0
    end do

! - Loop on Gauss points
    do ipg = 1, npg
        ldec = (ipg-1)*nbNode

! ----- Compute normal
        nx = 0.d0
        ny = 0.d0
        call vff2dn(cellDime, nbNode, ipg, jvWeight, jvDShape,&
                    zr(jvGeom), nx, ny, poids)
        if (l_axis) then
            r = 0.d0
            do i = 1, nbNode
                r = r + zr(jvGeom+2*(i-1))*zr(jvShape+ldec+i-1)
            end do
            poids = poids*r
        endif

! ----- Get value of speed
        call evalFaceSpeedVale(lFunc    , lTime   , time  ,&
                               nbNode   , cellDime, ipg   ,&
                               jvShape  , jvGeom  , jvLoad,&
                               speedVale, x, y)

! ----- Get direction of speed
        call evalFaceSpeedDire(cellDime , jvLoad, speedDire, &
                               ipg, nx, ny,&
                               lFunc_ = lFunc, lReal_ = lReal,&
                               lTime_ = lTime, time_ = time ,&
                               x_ = x, y_ = y)

! ----- Compute vector
        do i = 1, nbNode
            ii = ndofbynode*i
            zr(jvVect+ii-1) = zr(jvVect+ii-1) - speedDire*&
                               poids *&
                               zr(jvShape+ldec+i-1) * speedVale * rho
        end do

    end do
!
end subroutine
