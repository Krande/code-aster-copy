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
subroutine te0179(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/vff2dn.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: ACOU / PLAN (boundary)
!
! Options: CHAR_ACOU_VFAC
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvGeom, jvMate, jvLoad, jvVect
    real(kind=8) :: rho, poids, nx, ny
    complex(kind=8) :: speedVale
    integer(kind=8) :: jvWeight, jvShape, jvDShape
    integer(kind=8) :: nbNode, npg, cellDime, ndof
    integer(kind=8) :: ldec
    integer(kind=8) :: i, ipg
    aster_logical :: l_axis
    real(kind=8) :: r
    integer(kind=8) :: j_mater
!
! --------------------------------------------------------------------------------------------------
!

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMate)
    call jevech('PVITEFC', 'L', jvLoad)

! - Get element parameters
    l_axis = (lteatt('AXIS', 'OUI'))
    call elrefe_info(fami='RIGI', &
                     nno=nbNode, npg=npg, ndim=cellDime, &
                     jpoids=jvWeight, jvf=jvShape, jdfde=jvDShape)
    ndof = nbNode

! - Get material properties
    j_mater = zi(jvMate)
    call getFluidPara(j_mater, rho_=rho)

! - Output field
    call jevech('PVECTTC', 'E', jvVect)
    do i = 1, ndof
        zc(jvVect+i-1) = (0.d0, 0.d0)
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
        speedVale = zc(jvLoad-1+1)

! ----- Compute vector
        do i = 1, nbNode
            zc(jvVect+i-1) = zc(jvVect+i-1)+ &
                             poids* &
                             zr(jvShape+ldec+i-1)*speedVale*rho
        end do

    end do
!
end subroutine
