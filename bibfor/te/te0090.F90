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
subroutine te0090(option, nomte)
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/vff2dn.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tefrep.h"
!
character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MECA_2D (skin)
!
! Options: CHAR_MECA_FR1D2D
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: ndimSpace = 2, nDimCell = 1
    integer :: jvWeight, jvShape, jvDShape
    integer :: jvGeom, jvForc, jvVect
    integer :: nno, npg
    integer :: kpg, iNode, iDof
    integer :: jdec, kdec
    real(kind=8) :: jacWeight, r, fx, fy, nx, ny
    aster_logical :: lAxis
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI',&
                     nno=nno, npg=npg,&
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShape)
    lAxis = lteatt('AXIS','OUI')

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call tefrep(option, 'PFR1D2D', jvForc)

! - Get output fields
    call jevech('PVECTUR', 'E', jvVect)
    do iDof = 1, ndimSpace*nno
        zr(jvVect+iDof-1) = 0.d0
    end do

! - Loop on Gauss points
    do kpg = 1, npg
        kdec = (kpg-1)*nno

! ----- Compute jacobian
        call vff2dn(nDimCell, nno, kpg, jvWeight, jvDShape,&
                    zr(jvGeom), nx, ny, jacWeight)

        if (lAxis) then
            r = 0.d0
            do iNode = 1, nno
                r = r + zr(jvGeom+ndimSpace*(iNode-1))*zr(jvShape+kdec+iNode-1)
            end do
            jacWeight = jacWeight * r
        endif

! ----- Compute force at Gauss point from node value
        fx = 0.d0
        fy = 0.d0
        do iNode = 1, nno
            jdec = (iNode-1) * ndimSpace
            fx = fx + zr(jvShape+kdec+iNode-1) * zr(jvForc+jdec )
            fy = fy + zr(jvShape+kdec+iNode-1) * zr(jvForc+jdec+1)
        end do

! ----- Compute force 
        do iNode = 1, nno
            zr(jvVect+ndimSpace*(iNode-1))   = zr(jvVect+ndimSpace*(iNode-1)) +&
                                               jacWeight * fx * zr(jvShape+kdec+iNode-1)
            zr(jvVect+ndimSpace*(iNode-1)+1) = zr(jvVect+ndimSpace*(iNode-1)+1) +&
                                               jacWeight * fy * zr(jvShape+kdec+iNode-1)
        end do
    end do
!
end subroutine
