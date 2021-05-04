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
subroutine te0173(option, nomte)
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
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_FLUIDE (boundary)
!
! Options: CHAR_MECA_VFAC
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lFunc, lTime
    integer :: jvGeom, jvMate, jvLoad, jvTime, jvVect
    aster_logical :: lReal
    real(kind=8) :: x, y, z, speedDire
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: jac, rho
    real(kind=8) :: time, speedVale
    integer :: jvWeight, jvShape, jvDShapeX, jvDShapeY
    integer :: nbNode, npg, cellDime, ndofbynode
    integer :: idec, jdec, kdec, ldec
    integer :: i, ii, ino, j, jno, ipg
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
    time  = 0.d0
    if (jvTime .ne. 0) then
        lTime = ASTER_TRUE
        time  = zr(jvTime)
    endif

! - Get element parameters
    call teattr('S', 'FORMULATION', fsi_form, iret)
    call elrefe_info(fami='RIGI',&
                     nno=nbNode, npg=npg, ndim=cellDime,&
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShapeX)
    ASSERT(nbNode .le. 9)
    jvDShapeY = jvDShapeX + 1
    if (fsi_form .eq. 'FSI_UPPHI') then
        ndofbynode = 2
    elseif (fsi_form .eq. 'FSI_UP' .or. fsi_form .eq. 'FSI_UPSI') then
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

! - CALCUL DES PRODUITS VECTORIELS OMI X OMJ
    do ino = 1, nbNode
        i = jvGeom + 3*(ino-1) -1
        do jno = 1, nbNode
            j = jvGeom + 3*(jno-1) -1
            sx(ino,jno) = zr(i+2) * zr(j+3) - zr(i+3) * zr(j+2)
            sy(ino,jno) = zr(i+3) * zr(j+1) - zr(i+1) * zr(j+3)
            sz(ino,jno) = zr(i+1) * zr(j+2) - zr(i+2) * zr(j+1)
        end do
    end do

! - Loop on Gauss points
    do ipg = 1, npg
        kdec = (ipg-1)*nbNode*cellDime
        ldec = (ipg-1)*nbNode

! ----- Compute normal
        nx = 0.d0
        ny = 0.d0
        nz = 0.d0
        do i = 1, nbNode
            idec = (i-1)*cellDime
            do j = 1, nbNode
                jdec = (j-1)*cellDime
                nx = nx + zr(jvDShapeX+kdec+idec) * zr(jvDShapeY+kdec+jdec) * sx(i,j)
                ny = ny + zr(jvDShapeX+kdec+idec) * zr(jvDShapeY+kdec+jdec) * sy(i,j)
                nz = nz + zr(jvDShapeX+kdec+idec) * zr(jvDShapeY+kdec+jdec) * sz(i,j)
            end do
        end do

! ----- Compute jacobian
        jac = sqrt (nx*nx + ny*ny + nz*nz)

! ----- Get value of speed
        call evalFaceSpeedVale(lFunc    , lTime   , time  ,&
                               nbNode   , cellDime, ipg   ,&
                               jvShape  , jvGeom  , jvLoad,&
                               speedVale, x, y, z)

! ----- Get direction of speed
        call evalFaceSpeedDire(fsi_form, cellDime , jvLoad, speedDire, &
                               ipg      , nx    , ny      ,&
                               lFunc_ = lFunc, lReal_ = lReal,&
                               lTime_ = lTime, time_ = time ,&
                               x_ = x, y_ = y,&
                               z_ = z, nz_ = nz)

! ----- Compute vector
        do i = 1, nbNode
            ii = ndofbynode*i
            zr(jvVect+ii-1) = zr(jvVect+ii-1) - speedDire *&
                               jac*zr(jvWeight+ipg-1) *&
                               zr(jvShape+ldec+i-1) * speedVale * rho
        end do

    end do
!
end subroutine
