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
subroutine te0182(option, nomte)
!.......................................................................
!
!*    BUT: CALCUL DES MATRICES ELEMENTAIRES EN ACOUSTIQUE
!*         CORRESPONDANT AU TERME D'AMORTISSEMENT ACOUSTIQUE
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!*         OPTION : 'AMOR_ACOU'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    implicit none
!*
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/getFluidPara.h"
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: nomte, option
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), jac
    complex(kind=8) :: rhosz, cele_c
    real(kind=8) :: rho, alpha, q_alpha, cele_r, cele_i, onde_flui
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom, imate, jv_amor
    integer(kind=8) :: ndim, nno, ndi, ipg, imattt
    integer(kind=8) :: idec, jdec, kdec, ldec, nnos, npg2, jgano
    integer(kind=8) :: i, ij, ino, j, jno, mater
!
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
    idfdy = idfdx+1
    ndi = nno*(nno+1)/2
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATTTC', 'E', imattt)
    call jevech('PMATERC', 'L', imate)
    call jevech('PWATFLAC', 'L', jv_amor)
!
! - Get material properties for fluid
!
    mater = zi(imate)

    call getFluidPara(mater, rho_=rho, cele_r_=cele_r, cele_i_=cele_i, alpha_=alpha)

    cele_c = dcmplx(cele_r, cele_i)
!
! - Conditions on fluid parameters
!
    if ((1.d0-alpha) .ge. r8prem()) then
        q_alpha = (1.d0+alpha)/(1.d0-alpha)
    else
        call utmess('F', 'FLUID1_5', sk='ALPHA')
    end if

    if (rho .le. r8prem()) then
        call utmess('F', 'FLUID1_6', sk='RHO')
    end if

    if ((abs(cele_r) .le. r8prem()) .and. (abs(cele_i) .le. r8prem())) then
        call utmess('F', 'FLUID1_7', sk='CELE_R + i*CELE_I')
    end if

    if (zi(jv_amor-1+1) .eq. 1) then
        onde_flui = +1.d0
    else
        onde_flui = -1.d0
    end if
!
! - Output field
!

    rhosz = 1.d0/(q_alpha*cele_c)

!
! - Compute
!
    do i = 1, ndi
        zc(imattt+i-1) = (0.0d0, 0.0d0)
    end do
!
!    CALCUL DES PRODUITS VECTORIELS OMI X OMJ
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg2
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
!
!   CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
!
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
!
                nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
!
            end do
        end do
!
!   CALCUL DU JACOBIEN AU POINT DE GAUSS IPG
!
        jac = sqrt(nx*nx+ny*ny+nz*nz)
!
        do i = 1, nno
            do j = 1, i
                ij = (i-1)*i/2+j
                zc(imattt+ij-1) = zc(imattt+ij-1)+jac*onde_flui*rhosz*&
                                          &zr(ipoids+ipg-1)*zr(ivf+ldec+i-&
                                          &1)*zr(ivf+ldec+j-1)
            end do
        end do
    end do
end subroutine
