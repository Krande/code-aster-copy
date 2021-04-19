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
subroutine te0028(option, nomte)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tefrep.h"
!
character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MECA_3D (skin)
!
! Options: CHAR_MECA_FR2D3D
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: ndimSpace = 3, ndimCell = 2
    integer :: jvWeight, jvShape, jvDShapeX, jvDShapeY
    integer :: jvGeom, jvForc, jvVect
    integer :: nno, npg
    integer :: kpg, iNode, iDof, jNode
    integer :: kdec, ldec, idec, jdec
    real(kind=8) :: jacWeight, fx, fy, fz
    real(kind=8) :: jac, nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI',&
                     nno=nno, npg=npg,&
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShapeX)
    jvDShapeY = jvDShapeX + 1

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call tefrep(option, 'PFR2D3D', jvForc)

! - Get output fields
    call jevech('PVECTUR', 'E', jvVect)
    do iDof = 1, ndimSpace*nno
        zr(jvVect+iDof-1) = 0.d0
    end do

! - Compute exterior product of face
    do iNode = 1, nno
        idec = jvGeom + 3*(iNode-1) -1
        do jNode = 1, nno
            jdec = jvGeom + 3*(jNode-1) -1
            sx(iNode,jNode) = zr(idec+2) * zr(jdec+3) - zr(idec+3) * zr(jdec+2)
            sy(iNode,jNode) = zr(idec+3) * zr(jdec+1) - zr(idec+1) * zr(jdec+3)
            sz(iNode,jNode) = zr(idec+1) * zr(jdec+2) - zr(idec+2) * zr(jdec+1)
        end do
    end do

! - Loop on Gauss points
    do kpg = 1, npg
        kdec = (kpg-1)*nno*ndimCell
        ldec = (kpg-1)*nno

! ----- Compute normal to face
        nx = 0.d0
        ny = 0.d0
        nz = 0.d0
        do iNode = 1, nno
            idec = (iNode-1)*ndimCell
            do jNode = 1, nno
                jdec = (jNode-1)*ndimCell
                nx = nx + zr(jvDShapeX+kdec+idec) * zr(jvDShapeY+kdec+jdec) * sx(iNode,jNode)
                ny = ny + zr(jvDShapeX+kdec+idec) * zr(jvDShapeY+kdec+jdec) * sy(iNode,jNode)
                nz = nz + zr(jvDShapeX+kdec+idec) * zr(jvDShapeY+kdec+jdec) * sz(iNode,jNode)
            end do
        end do

! ----- Compute jacobian
        jac       = sqrt(nx*nx + ny*ny + nz*nz)
        jacWeight = jac * zr(jvWeight+kpg-1)

! ----- Compute force at Gauss point from node value
        fx = 0.d0
        fy = 0.d0
        fz = 0.d0
        do iNode = 1, nno
            kdec = ndimSpace*(iNode-1)
            fx = fx + zr(jvShape+ldec+iNode-1) * zr(jvForc+kdec)
            fy = fy + zr(jvShape+ldec+iNode-1) * zr(jvForc+kdec+1)
            fz = fz + zr(jvShape+ldec+iNode-1) * zr(jvForc+kdec+2)
        end do
        do iNode = 1, nno
            kdec = ndimSpace*(iNode-1)
            zr(jvVect+kdec  ) = zr(jvVect+kdec  ) + jacWeight * fx * zr(jvShape+ldec+iNode-1)
            zr(jvVect+kdec+1) = zr(jvVect+kdec+1) + jacWeight * fy * zr(jvShape+ldec+iNode-1)
            zr(jvVect+kdec+2) = zr(jvVect+kdec+2) + jacWeight * fz * zr(jvShape+ldec+iNode-1)
        end do
    end do
!
end subroutine
