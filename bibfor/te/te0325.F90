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
subroutine te0325(option, nomte)
    implicit none
!
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES DE FLUX FLUIDE EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 2D
!
!          OPTION : 'CHAR_THER_ACCE_R 'OU 'CHAR_THER_ACCE_X'
!                    OU 'CHAR_THER_ACCE_Y'OU 'CHAR_THER_ACCE_Z'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
!
    integer(kind=8) :: icodre(1), spt
    character(len=8) :: fami, poum
    character(len=16) :: nomte, option
    real(kind=8) :: jac, nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: norm(3)
    real(kind=8) :: acloc(3, 9), acc(3, 9), flufn(9)
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom
    integer(kind=8) :: ndim, nno, ipg, npg1, ivectt, imate
    integer(kind=8) :: idec, jdec, kdec, ldec, nnos, jgano
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacce, idim, ino, itemp, j, jno
    integer(kind=8) :: k, mater
    real(kind=8) :: rho(1)
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PVECTTR', 'E', ivectt)
!
    fami = 'RIGI'
    spt = 1
    poum = '+'
    mater = zi(imate)
!
    if (option(16:16) .eq. 'R') then
        call jevech('PACCELR', 'L', iacce)
    elseif (option(16:16) .eq. 'X' .or. option(16:16) .eq. 'Y' .or. &
            option(16:16) .eq. 'Z') then
        call jevech('PTEMPER', 'L', itemp)
    end if
!
! ON RECUPERE LE CHAMNO DE DEPL (MODAL)
!
    k = 0
    do i = 1, nno
        if (option(16:16) .eq. 'R') then
            do idim = 1, 3
                k = k+1
                acloc(idim, i) = zr(iacce+k-1)
            end do
        else if (option(16:16) .eq. 'X') then
            k = k+1
            acloc(1, i) = zr(itemp+k-1)
            acloc(2, i) = 0.d0
            acloc(3, i) = 0.d0
        else if (option(16:16) .eq. 'Y') then
            k = k+1
            acloc(1, i) = 0.d0
            acloc(2, i) = zr(itemp+k-1)
            acloc(3, i) = 0.d0
        else if (option(16:16) .eq. 'Z') then
            k = k+1
            acloc(1, i) = 0.d0
            acloc(2, i) = 0.d0
            acloc(3, i) = zr(itemp+k-1)
        end if
    end do
!
    do i = 1, nno
        zr(ivectt+i-1) = 0.d0
    end do
!
!     CALCUL DES PRODUITS VECTORIELS OMI X OMJ
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
!     BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        call rcvalb(fami, ipg, spt, poum, mater, &
                    ' ', 'THER', 0, ' ', [0.d0], &
                    1, 'RHO_CP', rho, icodre, 1)
!
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
!        --- ON CALCULE L ACCEL AU POINT DE GAUSS
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
        acc(1, ipg) = 0.0d0
        acc(2, ipg) = 0.0d0
        acc(3, ipg) = 0.0d0
        do i = 1, nno
            acc(1, ipg) = acc(1, ipg)+acloc(1, i)*zr(ivf+ldec+i-1)
            acc(2, ipg) = acc(2, ipg)+acloc(2, i)*zr(ivf+ldec+i-1)
            acc(3, ipg) = acc(3, ipg)+acloc(3, i)*zr(ivf+ldec+i-1)
        end do
!
!        CALCUL DU JACOBIEN AU POINT DE GAUSS IPG
!
        jac = sqrt(nx*nx+ny*ny+nz*nz)
        norm(1) = nx/jac
        norm(2) = ny/jac
        norm(3) = nz/jac
        flufn(ipg) = 0.d0
! CALCUL DU FLUX FLUIDE NORMAL AU POINT DE GAUSS
        flufn(ipg) = acc(1, ipg)*norm(1)+acc(2, ipg)*norm(2)+acc(3, ipg)*norm(3)
!
        do i = 1, nno
            zr(ivectt+i-1) = zr(ivectt+i-1)+jac*zr(ipoids+ipg-1)*flufn(ipg)*rho(1)*zr(ivf+&
                             &ldec+i-1)
        end do
    end do
!
end subroutine
