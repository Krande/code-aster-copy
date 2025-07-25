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
subroutine shl329()
    implicit none
!....................................................................
!   CALCUL DES TERMES ELEMENTAIRES DE L'ACCEPTANCE
!     OPTION : ACCEPTANCE
!....................................................................
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/wkvect.h"
    character(len=7) :: ielem, imode
    character(len=24) :: vetel
    real(kind=8) :: sx(9, 9), sy(9, 9), sz(9, 9), jac(9)
    real(kind=8) :: nx(9), ny(9), nz(9), norm(3, 9), acc(3, 9)
    real(kind=8) :: flufn(9), acloc(3, 8), qsi, eta, zero, un
    real(kind=8) :: x(3, 9), ff(4, 4), dfdx(4, 4), dfdy(4, 4)
    integer(kind=8) :: igeom, ipg, iadzi, iazk24, ivetel
    integer(kind=8) :: iacce, iharm, ivectu, i, k, idim, ino, jno, j
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
!
! DEB ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
    call jevech('PACCELR', 'L', iacce)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PNUMMOD', 'L', iharm)
    call jevech('PVECTUR', 'E', ivectu)
!
    zero = 0.0d0
    un = 1.0d0
!
    do i = 1, nno
        acloc(1, i) = zero
        acloc(2, i) = zero
        acloc(3, i) = zero
    end do
!
    k = 0
    do i = 1, nno
        do idim = 1, 3
            k = k+1
            acloc(idim, i) = zr(iacce+k-1)
        end do
    end do
!
    do ipg = 1, npg
        acc(1, ipg) = zero
        acc(2, ipg) = zero
        acc(3, ipg) = zero
    end do
!
    do ipg = 1, npg
!
        qsi = zr(icoopg-1+ndim*(ipg-1)+1)
        eta = zr(icoopg-1+ndim*(ipg-1)+2)
!
        ff(1, ipg) = (un-qsi)*(un-eta)/4.d0
        ff(2, ipg) = (un+qsi)*(un-eta)/4.d0
        ff(3, ipg) = (un+qsi)*(un+eta)/4.d0
        ff(4, ipg) = (un-qsi)*(un+eta)/4.d0
!
        dfdx(1, ipg) = -(un-eta)/4.d0
        dfdx(2, ipg) = (un-eta)/4.d0
        dfdx(3, ipg) = (un+eta)/4.d0
        dfdx(4, ipg) = -(un+eta)/4.d0
!
        dfdy(1, ipg) = -(un-qsi)/4.d0
        dfdy(2, ipg) = -(un+qsi)/4.d0
        dfdy(3, ipg) = (un+qsi)/4.d0
        dfdy(4, ipg) = (un-qsi)/4.d0
!
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
    do ipg = 1, npg
        nx(ipg) = zero
        ny(ipg) = zero
        nz(ipg) = zero
        do i = 1, nno
            do j = 1, nno
                nx(ipg) = nx(ipg)+dfdx(i, ipg)*dfdy(j, ipg)*sx(i, j)
                ny(ipg) = ny(ipg)+dfdx(i, ipg)*dfdy(j, ipg)*sy(i, j)
                nz(ipg) = nz(ipg)+dfdx(i, ipg)*dfdy(j, ipg)*sz(i, j)
            end do
        end do
!
!      CALCUL DU JACOBIEN AU POINT DE GAUSS IPG
!
        jac(ipg) = sqrt(nx(ipg)*nx(ipg)+ny(ipg)*ny(ipg)+nz(ipg)*nz(ipg))
!
!       CALCUL DE LA NORMALE UNITAIRE
!
        norm(1, ipg) = nx(ipg)/jac(ipg)
        norm(2, ipg) = ny(ipg)/jac(ipg)
        norm(3, ipg) = nz(ipg)/jac(ipg)
    end do
!
!
    do ipg = 1, npg
        do i = 1, nno
            acc(1, ipg) = acc(1, ipg)+acloc(1, i)*ff(i, ipg)
            acc(2, ipg) = acc(2, ipg)+acloc(2, i)*ff(i, ipg)
            acc(3, ipg) = acc(3, ipg)+acloc(3, i)*ff(i, ipg)
        end do
    end do
!
!    CALCUL DE COORDONNEES AUX POINTS DE GAUSS
!
    do ipg = 1, npg
        x(1, ipg) = zero
        x(2, ipg) = zero
        x(3, ipg) = zero
        do j = 1, nno
            x(1, ipg) = x(1, ipg)+zr(igeom+3*(j-1)-1+1)*ff(j, ipg)
            x(2, ipg) = x(2, ipg)+zr(igeom+3*(j-1)-1+2)*ff(j, ipg)
            x(3, ipg) = x(3, ipg)+zr(igeom+3*(j-1)-1+3)*ff(j, ipg)
        end do
!
! CALCUL DU FLUX FLUIDE NORMAL AUX POINTS DE GAUSS
!
        flufn(ipg) = acc(1, ipg)*norm(1, ipg)+acc(2, ipg)*norm(2, ipg)+acc(3, ipg)*norm(3, ipg)
!
    end do
!
! STOCKAGE DU FLUX FLUIDE DANS UN VECTEUR INDEXE
! PAR LE MODE ET L'ELEMENT
!
    imode = 'CHBIDON'
    ielem = 'CHBIDON'
    call codent(zi(iharm), 'D0', imode)
    call tecael(iadzi, iazk24, noms=0)
    call codent(zi(iadzi), 'D0', ielem)
    vetel = '&&329.M'//imode//'.EL'//ielem
!        ON CONSERVE L'ALLOCATION DYNAMIQUE AU DETRIMENT DE L'ALLOCATION
!        STATIQUE, CAR VETEL EST UTILIE A L'EXTERIEUR DES ROUTINES
!        ELEMENTAIRES
    call wkvect(vetel, 'V V R8', 4*npg, ivetel)
    do ipg = 0, npg-1
        zr(ivetel+4*ipg) = jac(ipg+1)*flufn(ipg+1)
        zr(ivetel+4*ipg+1) = x(1, ipg+1)
        zr(ivetel+4*ipg+2) = x(2, ipg+1)
        zr(ivetel+4*ipg+3) = x(3, ipg+1)
    end do
!
end subroutine
