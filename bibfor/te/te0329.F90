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
subroutine te0329(option, nomte)
    implicit none
!....................................................................
!   CALCUL DES TERMES ELEMENTAIRES DE L'ACCEPTANCE
!     OPTION : ACCEPTANCE
!....................................................................
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/shl329.h"
#include "asterfort/tecael.h"
#include "asterfort/wkvect.h"
    character(len=7) :: ielem, imode
    character(len=16) :: nomte, option
    character(len=24) :: vetel
    real(kind=8) :: sx(9, 9), sy(9, 9), sz(9, 9), jac(9)
    real(kind=8) :: nx(9), ny(9), nz(9), norm(3, 9), acc(3, 9)
    real(kind=8) :: flufn(9), acloc(3, 8)
    real(kind=8) :: x(3, 9)
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom
    integer(kind=8) :: ndim, nno, ipg, npg1, npg4
    integer(kind=8) :: idec, jdec, kdec, ldec
    integer(kind=8) :: nnos, jgano
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacce, iadzi, iazk24, idim, iharm, ino
    integer(kind=8) :: ivectu, ivetel, j, jno, k
!-----------------------------------------------------------------------
    if (lteatt('DIM_TOPO_MODELI', '3')) then
!   ----------------------------------------
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
        idfdy = idfdx+1
!
        ASSERT((npg1 == 3) .or. (npg1 == 4))
        call jevech('PACCELR', 'L', iacce)
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PNUMMOD', 'L', iharm)
        call jevech('PVECTUR', 'E', ivectu)
!
        do i = 1, nno
            acloc(1, i) = 0.0d0
            acloc(2, i) = 0.0d0
            acloc(3, i) = 0.0d0
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
        do ipg = 1, npg1
            acc(1, ipg) = 0.0d0
            acc(2, ipg) = 0.0d0
            acc(3, ipg) = 0.0d0
        end do
!
!
        do ipg = 1, npg1
            ldec = (ipg-1)*nno
            do i = 1, nno
                acc(1, ipg) = acc(1, ipg)+acloc(1, i)*zr(ivf+ldec+i-1)
                acc(2, ipg) = acc(2, ipg)+acloc(2, i)*zr(ivf+ldec+i-1)
                acc(3, ipg) = acc(3, ipg)+acloc(3, i)*zr(ivf+ldec+i-1)
            end do
        end do
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
!
!     BOUCLE SUR LES POINTS DE GAUSS
!
        do ipg = 1, npg1
!
            kdec = (ipg-1)*nno*ndim
            ldec = (ipg-1)*nno
!
            nx(ipg) = 0.0d0
            ny(ipg) = 0.0d0
            nz(ipg) = 0.0d0
!
            do i = 1, nno
                idec = (i-1)*ndim
                do j = 1, nno
                    jdec = (j-1)*ndim
!
                    nx(ipg) = nx(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                    ny(ipg) = ny(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                    nz(ipg) = nz(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
!
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
!    CALCUL DE COORDONNEES AUX POINTS DE GAUSS
!
        do ipg = 1, npg1
            ldec = (ipg-1)*nno
            x(1, ipg) = 0.0d0
            x(2, ipg) = 0.0d0
            x(3, ipg) = 0.0d0
!
            do j = 1, nno
!
                x(1, ipg) = x(1, ipg)+zr(igeom+3*(j-1)-1+1)*zr(ivf+ &
                                                               ldec+j-1)
                x(2, ipg) = x(2, ipg)+zr(igeom+3*(j-1)-1+2)*zr(ivf+ &
                                                               ldec+j-1)
                x(3, ipg) = x(3, ipg)+zr(igeom+3*(j-1)-1+3)*zr(ivf+ &
                                                               ldec+j-1)
!
            end do
!
! CALCUL DU FLUX FLUIDE NORMAL AUX POINTS DE GAUSS
!
            flufn(ipg) = acc(1, ipg)*norm(1, ipg)+acc(2, ipg)*norm(2, ipg)+acc(3, ipg)*norm(3, ipg)
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
!       on écrit toujours le résultat comme si l'élément avait 4 points
!       de Gauss
        npg4 = 4
        call wkvect(vetel, 'V V R8', 4*npg4, ivetel)
        do ipg = 0, npg1-1
            zr(ivetel+4*ipg) = jac(ipg+1)*zr(ipoids+ipg)*flufn(ipg+1)
            zr(ivetel+4*ipg+1) = x(1, ipg+1)
            zr(ivetel+4*ipg+2) = x(2, ipg+1)
            zr(ivetel+4*ipg+3) = x(3, ipg+1)
        end do
        if (npg1 == 3) then
            ipg = 3
            zr(ivetel+4*ipg) = 0.D0
            zr(ivetel+4*ipg+1) = 0.D0
            zr(ivetel+4*ipg+2) = 0.D0
            zr(ivetel+4*ipg+3) = 0.D0
        end if
!
!
    else if (nomte .eq. 'MEDKQU4') then
!   ----------------------------------
        call shl329()
!
    end if
!
end subroutine
