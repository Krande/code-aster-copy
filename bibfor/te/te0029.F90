! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine te0029(option, nomte)
    implicit none
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          CORRESPONDANT A UN CHARGEMENT FORCE_FACE FONCTION
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_MECA_FF2D3D '
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    character(len=8) :: nompar(4)
    character(len=16) :: nomte, option
    real(kind=8) :: jac, nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: valpar(4)
    integer :: ipoids, ivf, idfdx, idfdy, igeom, itemps
    integer :: ndim, nno, ipg, npg1, iforc, nnos
    integer :: idec, jdec, kdec, ldec, jgano
!
!-----------------------------------------------------------------------
    integer :: i, ier, ii, ino, ires, j, jno
!
    real(kind=8) :: fx, fy, fz, xx, yy, zz
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx + 1
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PFF2D3D', 'L', iforc)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PVECTUR', 'E', ires)
!
    valpar(4) = zr(itemps)
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
    nompar(4) = 'INST'
!
    do i = 1, 3*nno
        zr(ires+i-1) = 0.0d0
    end do
!
!    CALCUL DES PRODUITS VECTORIELS OMI X OMJ
!
    do ino = 1, nno
        i = igeom + 3*(ino-1) -1
        do jno = 1, nno
            j = igeom + 3*(jno-1) -1
            sx(ino,jno) = zr(i+2) * zr(j+3) - zr(i+3) * zr(j+2)
            sy(ino,jno) = zr(i+3) * zr(j+1) - zr(i+1) * zr(j+3)
            sz(ino,jno) = zr(i+1) * zr(j+2) - zr(i+2) * zr(j+1)
        end do
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
!    CALCUL DE FX,FY,FZ
!
        xx = 0.d0
        yy = 0.d0
        zz = 0.d0
        do i = 1, nno
            xx = xx + zr(igeom+3*(i-1)) * zr(ivf+ldec+i-1)
            yy = yy + zr(igeom+3*i-2) * zr(ivf+ldec+i-1)
            zz = zz + zr(igeom+3*i-1) * zr(ivf+ldec+i-1)
        end do
        valpar(1) = xx
        valpar(2) = yy
        valpar(3) = zz
        call fointe('FM', zk8(iforc), 4, nompar, valpar,&
                    fx, ier)
        call fointe('FM', zk8(iforc+1), 4, nompar, valpar,&
                    fy, ier)
        call fointe('FM', zk8(iforc+2), 4, nompar, valpar,&
                    fz, ier)
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
                nx = nx + zr(idfdx+kdec+idec) * zr(idfdy+kdec+jdec) * sx(i,j)
                ny = ny + zr(idfdx+kdec+idec) * zr(idfdy+kdec+jdec) * sy(i,j)
                nz = nz + zr(idfdx+kdec+idec) * zr(idfdy+kdec+jdec) * sz(i,j)
!
            end do
        end do
!
!   LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE
!
        jac = sqrt (nx*nx + ny*ny + nz*nz)
!
        do i = 1, nno
            ii = 3 * (i-1) -1
!
            zr(ires+ii+1) = zr(ires+ii+1) + zr(ipoids+ipg-1) * fx * zr(ivf+ldec+i-1) * jac
!
            zr(ires+ii+2) = zr(ires+ii+2) + zr(ipoids+ipg-1) * fy * zr(ivf+ldec+i-1) * jac
!
            zr(ires+ii+3) = zr(ires+ii+3) + zr(ipoids+ipg-1) * fz * zr(ivf+ldec+i-1) * jac
!
        end do
!
    end do
end subroutine
