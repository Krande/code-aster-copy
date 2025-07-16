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
subroutine te0374(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
#include "asterc/r8prem.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_FLUIDE (boundary)
!
! Options: ONDE_FLUI
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_geom, jv_mate, jv_onde, jv_matr
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: jac, rho, celer
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy
    integer(kind=8) :: nno, ndi, npg, ndim
    integer(kind=8) :: idec, jdec, kdec, ldec
    integer(kind=8) :: i, ii, ij, ino, j, jj, ipg
    integer(kind=8) :: jno
    integer(kind=8) :: j_mater, iret
    character(len=16) :: FEForm
!
! --------------------------------------------------------------------------------------------------
!

!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
    call jevech('PONDECR', 'L', jv_onde)
!
! - Get element parameters
!
    call teattr('S', 'FORMULATION', FEForm, iret)
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx)
    ASSERT(nno .le. 9)
    idfdy = idfdx+1
    if (FEForm .eq. 'U_P_PHI') then
        ndi = nno*(2*nno+1)
    elseif (FEForm .eq. 'U_P' .or. FEForm .eq. 'U_PSI') then
        ndi = nno*(nno+1)/2
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho_=rho, cele_r_=celer)
!
! - Output field
!
    call jevech('PMATUUR', 'E', jv_matr)
    do i = 1, ndi
        zr(jv_matr+i-1) = 0.d0
    end do
!
! - Compute
!
    if (abs(zr(jv_onde)) .gt. r8prem()) then
        do ino = 1, nno
            i = jv_geom+3*(ino-1)-1
            do jno = 1, nno
                j = jv_geom+3*(jno-1)-1
                sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
                sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
                sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
            end do
        end do
! ----- Loop on Gauss points
        do ipg = 1, npg
            kdec = (ipg-1)*nno*ndim
            ldec = (ipg-1)*nno
! --------- Compute normal
            nx = 0.0d0
            ny = 0.0d0
            nz = 0.0d0
            do i = 1, nno
                idec = (i-1)*ndim
                do j = 1, nno
                    jdec = (j-1)*ndim
                    nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                    ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                    nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
                end do
            end do
! --------- Compute jacobian
            jac = sqrt(nx*nx+ny*ny+nz*nz)
            if (FEForm .eq. 'U_P_PHI') then
                do i = 1, nno
                    do j = 1, i
                        ii = 2*i
                        jj = 2*j
                        ij = (ii-1)*ii/2+jj
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)- &
                                           jac*zr(ipoids+ipg-1)* &
                                           zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)*rho/celer
                    end do
                end do
            elseif (FEForm .eq. 'U_P') then
                do i = 1, nno
                    do j = 1, i
                        ij = (i-1)*i/2+j
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)+ &
                                           jac*zr(ipoids+ipg-1)* &
                                           zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)/celer
                    end do
                end do
            elseif (FEForm .eq. 'U_PSI') then
                do i = 1, nno
                    do j = 1, i
                        ij = (i-1)*i/2+j
                        zr(jv_matr+ij-1) = zr(jv_matr+ij-1)- &
                                           jac*zr(ipoids+ipg-1)* &
                                           zr(ivf+ldec+i-1)*zr(ivf+ldec+j-1)*rho/celer
                    end do
                end do
            else
                call utmess('F', 'FLUID1_2', sk=FEForm)
            end if
        end do
    end if
!
end subroutine
