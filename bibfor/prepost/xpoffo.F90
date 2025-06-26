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
subroutine xpoffo(ndim, ndime, elrefp, nnop, igeom, &
                  co, ff)
! aslint: disable=W1306
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
    integer(kind=8) :: ndim, ndime, nnop, igeom
    real(kind=8) :: co(ndim), ff(nnop)
    character(len=8) :: elrefp
!
! person_in_charge: samuel.geniaut at edf.fr
!
!     FF : FONCTIONS DE FORMES AU NOEUD SOMMET OU D'INTERSECTION
!
!   IN
!     NDIM   : DIMENSION DU MAILLAGE
!     NDIME  : DIMENSION TOPOLOGIQUE DE LA MAILLE PARENT
!     ELREFP : ÉLÉMENT DE RÉFÉRENCE PARENT
!     NNOP   : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
!     IGEOM  : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
!     CO     : COORDONNÉES DU NOUVEAU NOEUD
!
!   OUT
!     FF    : FONCTIONS DE FORMES DE L'ELEMENT PARENT AU NOUVEAU NOEUD
!
!
    real(kind=8) :: a(3), b(3), c(3), ab(3), ac(3), nd(3), norme, nab, y(3)
    real(kind=8) :: an(ndim), coloc(2), xe(3), n(ndim)
    integer(kind=8) :: j, igeolo, ino
    character(len=24) :: geomlo
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!
!     CAS DES ELEMENTS PRINCIPAUX : C SIMPLE, ON APPELLE REEREF
    if (ndim .eq. ndime) then
        call reeref(elrefp, nnop, zr(igeom), co, ndim, &
                    xe, ff)
!
!
!     CAS DES ELEMENTS DE BORDS : C PLUS COMPLIQUÉ
!     ON NE PROCEDE PAS COMME DANS TE0036 CAR ICI, ON N'EST PAS
!     DANS UNE BOUCLE SUR LES SOUS-ELEMENTS
    else if (ndim .ne. ndime) then
!
        ASSERT(ndim .eq. ndime+1)
!
!       CREATION D'UNE BASE LOCALE A LA MAILLE PARENT ABCD
        ab(:) = 0.d0
        ac(:) = 0.d0
        do j = 1, ndim
            a(j) = zr(igeom-1+ndim*(1-1)+j)
            b(j) = zr(igeom-1+ndim*(2-1)+j)
            if (ndim .eq. 3) c(j) = zr(igeom-1+ndim*(3-1)+j)
            ab(j) = b(j)-a(j)
            if (ndim .eq. 3) ac(j) = c(j)-a(j)
        end do
!
        if (ndime .eq. 2) then
!         CREATION DU REPERE LOCAL LO : (AB,Y)
            call provec(ab, ac, nd)
            call normev(nd, norme)
            call normev(ab, nab)
            call provec(nd, ab, y)
        else if (ndime .eq. 1) then
!         CREATION DU REPERE LOCAL 1D : AB/NAB
            call normev(ab, nab)
            nd(:) = 0.d0
            nd(1) = ab(2)
            nd(2) = -ab(1)
        end if
!
!       COORDONNÉES DES NOEUDS DE L'ELREFP DANS LE REPÈRE LOCAL
        geomlo = '&&TE0036.GEOMLO'
        call wkvect(geomlo, 'V V R', nnop*ndime, igeolo)
        do ino = 1, nnop
            do j = 1, ndim
                n(j) = zr(igeom-1+ndim*(ino-1)+j)
                an(j) = n(j)-a(j)
            end do
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(igeolo-1+ndime*(ino-1)+1) = ddot(b_n, an, b_incx, ab, b_incy)
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            if (ndime .eq. 2) zr(igeolo-1+ndime*(ino-1)+2) = ddot(b_n, an, b_incx, y, b_incy)
        end do
!
!       COORDONNÉES RÉELLES LOCALES DU POINT EN QUESTION : COLOC
        do j = 1, ndim
            an(j) = co(j)-a(j)
        end do
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        coloc(1) = ddot(b_n, an, b_incx, ab, b_incy)
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        if (ndime .eq. 2) coloc(2) = ddot(b_n, an, b_incx, y, b_incy)
        call reeref(elrefp, nnop, zr(igeolo), coloc, ndime, &
                    xe, ff)
!
        call jedetr(geomlo)
!
    end if
!
    call jedema()
!
end subroutine
