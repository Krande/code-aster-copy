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

subroutine te0278(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: nomte, option
!.......................................................................
!
!     BUT: CALCUL DES MATRICES ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT AU TERME D'ECHANGE ENTRE 2 PAROIS (FACE)
!          D'ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'RESI_THER_PARO_F'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!
!
    character(len=8) :: nompar(4)
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), xx, yy, zz
    real(kind=8) :: jac, tem, hechp, valpar(4)
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom
    integer(kind=8) :: ndim, nno, ipg, npg1, iveres, ihechp
    integer(kind=8) :: idec, jdec, kdec, ldec
    integer(kind=8) ::  nnos, jgano
!---- DEBUT-----------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ier, ino, itemp, itemps, j, jno
!
!-----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PHECHPF', 'L', ihechp)
    call jevech('PTEMPEI', 'L', itemp)
    call jevech('PRESIDU', 'E', iveres)
!
    valpar(4) = zr(itemps)
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
    nompar(4) = 'INST'
    do i = 1, 2*nno
        zr(iveres+i-1) = 0.0d0
    end do
!
!    CALCUL DES PRODUITS VECTORIELS OMI * OMJ
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
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
!    CALCUL DE HECHP
!
        xx = 0.d0
        yy = 0.d0
        zz = 0.d0
        do i = 1, nno
            xx = xx+zr(igeom+3*i-3)*zr(ivf+ldec+i-1)
            yy = yy+zr(igeom+3*i-2)*zr(ivf+ldec+i-1)
            zz = zz+zr(igeom+3*i-1)*zr(ivf+ldec+i-1)
        end do
        valpar(1) = xx
        valpar(2) = yy
        valpar(3) = zz
        call fointe('A', zk8(ihechp), 4, nompar, valpar, &
                    hechp, ier)
        ASSERT(ier .eq. 0)
!
!   CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
!
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
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
! --- CALCUL DU JACOBIEN AU POINT DE GAUSS IPG
!
        jac = sqrt(nx*nx+ny*ny+nz*nz)
        tem = 0.d0
        do i = 1, nno
            ldec = (ipg-1)*nno
            tem = tem+(zr(itemp+nno+i-1)-zr(itemp+i-1))*zr(ivf+ldec+i-1)
        end do
        do i = 1, nno
            zr(iveres+i-1) = zr(iveres+i-1)-jac*hechp*zr(ipoids+ipg-1)*zr(ivf+ldec+i-1) &
                            &*tem
            zr(iveres+nno+i-1) = zr(iveres+nno+i-1)+jac*hechp*zr(ipoids+ipg-1)*zr(ivf+lde&
                                 &c+i-1)*tem
        end do
    end do
end subroutine
