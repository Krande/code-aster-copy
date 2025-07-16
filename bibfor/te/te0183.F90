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
subroutine te0183(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterfort/getFluidPara.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: ACOU / 3D (boundary)
!
! Options: CHAR_ACOU_VFAC
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvGeom, jvMate, jvLoad, jvVect
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: rho, jac
    complex(kind=8) :: speedVale
    integer(kind=8) :: jvWeight, jvShape, jvDShapeX, jvDShapeY
    integer(kind=8) :: nbNode, npg, cellDime, ndof
    integer(kind=8) :: idec, jdec, kdec, ldec
    integer(kind=8) :: i, ino, j, jno, ipg
    integer(kind=8) :: j_mater
!
! --------------------------------------------------------------------------------------------------
!

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMate)
    call jevech('PVITEFC', 'L', jvLoad)

! - Get element parameters
    call elrefe_info(fami='RIGI', &
                     nno=nbNode, npg=npg, ndim=cellDime, &
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShapeX)
    ASSERT(nbNode .le. 9)
    jvDShapeY = jvDShapeX+1
    ndof = nbNode

! - Get material properties
    j_mater = zi(jvMate)
    call getFluidPara(j_mater, rho_=rho)

! - Output field
    call jevech('PVECTTC', 'E', jvVect)
    do i = 1, ndof
        zc(jvVect+i-1) = (0.d0, 0.d0)
    end do

! - CALCUL DES PRODUITS VECTORIELS OMI X OMJ
    do ino = 1, nbNode
        i = jvGeom+3*(ino-1)-1
        do jno = 1, nbNode
            j = jvGeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
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
                nx = nx+zr(jvDShapeX+kdec+idec)*zr(jvDShapeY+kdec+jdec)*sx(i, j)
                ny = ny+zr(jvDShapeX+kdec+idec)*zr(jvDShapeY+kdec+jdec)*sy(i, j)
                nz = nz+zr(jvDShapeX+kdec+idec)*zr(jvDShapeY+kdec+jdec)*sz(i, j)
            end do
        end do

! ----- Compute jacobian
        jac = sqrt(nx*nx+ny*ny+nz*nz)

! ----- Get value of speed
        speedVale = zc(jvLoad-1+1)

! ----- Compute vector
        do i = 1, nbNode
            zc(jvVect+i-1) = zc(jvVect+i-1)+ &
                             jac*zr(jvWeight+ipg-1)* &
                             zr(jvShape+ldec+i-1)*speedVale*rho
        end do

    end do
!
end subroutine
