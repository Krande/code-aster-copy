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
subroutine te0417(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/dxroep.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/vectan.h"
#include "asterfort/vectci.h"
!
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     OPTION : 'MASS_INER'              (ELEMENTS COQUE_3D)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: xi(3, 9), xg(3), ix2(3), ix1x2, ix1x3, ix2x3, matine(6)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    integer(kind=8) :: i, intsn, intsx, j, jgeom, k, l1
    integer(kind=8) :: l2, lcastr, lzi, lzr, nb1, nb2, npgsn
    real(kind=8) :: epais, epais2, epais3, rho, rnormc, volume, wgt
    real(kind=8) :: xx, xy, xz, yy, yz, zz
!
! --------------------------------------------------------------------------------------------------
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
    call jevech('PGEOMER', 'L', jgeom)
!
    do i = 1, nb2
        xi(1, i) = zr(jgeom+3*(i-1))
        xi(2, i) = zr(jgeom+3*(i-1)+1)
        xi(3, i) = zr(jgeom+3*(i-1)+2)
    end do
!
    call dxroep(rho, epais)
    epais2 = epais*epais
    epais3 = epais*epais2
!
    call jevech('PMASSINE', 'E', lcastr)
!
    call vectan(nb1, nb2, xi, zr(lzr), vecta, &
                vectn, vectpt)
!
    volume = 0.d0
!
    do k = 1, 3
        xg(k) = 0.d0
        ix2(k) = 0.d0
    end do
    ix1x2 = 0.d0
    ix1x3 = 0.d0
    ix2x3 = 0.d0
!
    do intsn = 1, npgsn
!
!     RNORMC EST LE DETERMINANT DE LA SURFACE MOYENNE
!
        call vectci(intsn, nb1, xi, zr(lzr), rnormc)
!
!     WGT= ZR(9-1+INTE) * ZR(LZR+126-1+INTSN)
!        =    1.D0      * ZR(LZR+126-1+INTSN)
        wgt = zr(lzr+126-1+intsn)
!
        volume = volume+epais*wgt*rnormc
!
!     CENTRE DE GRAVITE
!
        l1 = lzr-1+135
        intsx = 8*(intsn-1)
        l2 = l1+intsx
!
        wgt = wgt*rnormc
!
        do j = 1, nb1
            do k = 1, 3
                xg(k) = xg(k)+epais*wgt*zr(l2+j)*xi(k, j)
            end do
!
!     MOMENTS ET PRODUITS D'INERTIE
!
            do i = 1, nb1
                do k = 1, 3
                    ix2(k) = ix2(k)+epais*wgt*zr(l2+j)*xi(k, j)*zr(l2+i)* &
                             xi(k, i)+epais3/12.d0*wgt*zr(l2+j)*vectn(j, k)*zr( &
                             l2+i)*vectn(i, k)
                end do
!
                ix1x2 = ix1x2+epais*wgt*zr(l2+j)*xi(1, j)*zr(l2+i)*xi(2, i)+ &
                        epais3/12.d0*wgt*zr(l2+j)*vectn(j, 1)*zr(l2+i)* &
                        vectn(i, 2)
                ix1x3 = ix1x3+epais*wgt*zr(l2+j)*xi(1, j)*zr(l2+i)*xi(3, i)+ &
                        epais3/12.d0*wgt*zr(l2+j)*vectn(j, 1)*zr(l2+i)* &
                        vectn(i, 3)
                ix2x3 = ix2x3+epais*wgt*zr(l2+j)*xi(2, j)*zr(l2+i)*xi(3, i)+ &
                        epais3/12.d0*wgt*zr(l2+j)*vectn(j, 2)*zr(l2+i)* &
                        vectn(i, 3)
            end do
!
        end do
!
    end do
!
    matine(1) = rho*(ix2(2)+ix2(3))
    matine(2) = rho*ix1x2
    matine(3) = rho*(ix2(1)+ix2(3))
    matine(4) = rho*ix1x3
    matine(5) = rho*ix2x3
    matine(6) = rho*(ix2(1)+ix2(2))
!
    zr(lcastr) = rho*volume
    zr(lcastr+1) = xg(1)/volume
    zr(lcastr+2) = xg(2)/volume
    zr(lcastr+3) = xg(3)/volume
!
    xx = zr(lcastr+1)*zr(lcastr+1)
    yy = zr(lcastr+2)*zr(lcastr+2)
    zz = zr(lcastr+3)*zr(lcastr+3)
    xy = zr(lcastr+1)*zr(lcastr+2)
    xz = zr(lcastr+1)*zr(lcastr+3)
    yz = zr(lcastr+2)*zr(lcastr+3)
!
    zr(lcastr+4) = matine(1)-zr(lcastr)*(yy+zz)
    zr(lcastr+5) = matine(3)-zr(lcastr)*(xx+zz)
    zr(lcastr+6) = matine(6)-zr(lcastr)*(xx+yy)
    zr(lcastr+7) = matine(2)-zr(lcastr)*xy
    zr(lcastr+8) = matine(4)-zr(lcastr)*xz
    zr(lcastr+9) = matine(5)-zr(lcastr)*yz
!
end subroutine
