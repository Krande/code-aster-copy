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
! aslint: disable=W1306
!
subroutine xdelt2(elp, n, ndime, ksi, &
                  ptint, ndim, tabco, tabls, ipp, ip, &
                  delta)
!
    implicit none
!
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/matinv.h"
#include "asterfort/provec.h"
    integer(kind=8) :: ndime, ndim, ipp, ip, n(3)
    real(kind=8) :: ksi(ndim), delta(ndime), ptint(*), tabco(*), tabls(*)
    character(len=8) :: elp
!                 CALCUL DE LA QUANTITE A MINIMISER POUR LE CALCUL
!                    DES COORDONNEES DU PT MILIEU DE LA FISSURE
!
!     ENTREE
!       ELP     : TYPE DE L'ELEMENT
!       N       : LES INDICES DES NOEUX D'UNE FACE DANS L'ELEMENT PARENT
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       KSI     : COORDONNEES DE REFERENCE DU POINT
!       PTINT  : COORDONNÉES DES POINTS D'INTERSECTION
!       TABCO   : COORDONNEES DES NOEUDS DE L'ELEMENT
!       TABLS   : VALEUR DES LSN DES NOEUDS DE L'ELEMENT
!
!     SORTIE
!       DELTA   : QUANTITE A MINIMISER
!     ----------------------------------------------------------------
!
    integer(kind=8) :: nno
    real(kind=8) :: refcoo(3, MT_NNOMAX), ff(MT_NNOMAX), dff(3, MT_NNOMAX)
    integer(kind=8) :: i, j, k
    real(kind=8) :: p(ndim), m(ndim), nor(ndim)
    real(kind=8) :: pint1(ndim), pint2(ndim)
    real(kind=8) :: r(ndime), det, dx
    real(kind=8) :: jac(ndime, ndime), inv(ndime, ndime)
    real(kind=8) :: v1(3), v2(3), nf(3)
!
!
!......................................................................
!
! --- CALCUL DES FONCTIONS DE FORME ET DE LEUR DERIVEES EN UN POINT
! --- DANS LA MAILLE
!
!
    pint1(:) = 0.d0
    pint2(:) = 0.d0
!
    do i = 1, ndim
        pint1(i) = ptint(ndim*(ipp-1)+i)
        pint2(i) = ptint(ndim*(ip-1)+i)
    end do
!
!     CALCUL DES FONCTIONS DE FORME DE L'ELEMENT EN KSI
    call elrfvf(elp, ksi, ff, nno)
!
!     CALCUL DES DERIVEES FONCTIONS DE FORME DE L'ELEMENT EN KSI
    call elrfdf(elp, ksi, dff)
!
    r(:) = 0.d0
    jac(:, :) = 0.d0
! --- CALCUL DE R1,DERIVEES PREMIERES DE R1 EN KSI
! ---           R1 : LEVEL SET NORMALE
    do j = 1, nno
        r(1) = r(1)+ff(j)*tabls(j)
        do i = 1, ndime
            jac(1, i) = jac(1, i)+dff(i, j)*tabls(j)
        end do
    end do
!
    m(:) = 0.d0
    nor(:) = 0.d0
    p(:) = 0.d0
!
! --- CALCUL DE R2 EN KSI
! ---           R2 : PLAN MEDIATEUR PASSANT PAR LES 2 PTS INTER
!               SPECIFIQUE NDIM=2
!
!     CHAQUE POINT P(X,Y,Z) DU PLAN MEDIATEUR VERIFIE NOR*MP=0
!                   M=((XI+XII)/2,(YI+YII)/2,(ZI+ZII)/2)
!                   X = SUM(FF*TABCO(1))
!                   Y = SUM(FF*TABCO(1))
!                   Z = SUM(FF*TABCO(1))
!                   NOR = (XII-XI,YII-YI,ZII-ZI)
!
    do i = 1, ndim
        do j = 1, nno
            p(i) = p(i)+ff(j)*tabco(ndim*(j-1)+i)
        end do
        m(i) = 0.5d0*(pint1(i)+pint2(i))
        nor(i) = pint2(i)-pint1(i)
        r(2) = r(2)+(p(i)-m(i))*nor(i)
    end do
!
!     CALCUL DES DERIVEES PREMIERES DE R2 EN KSI
    do i = 1, ndime
        do k = 1, ndim
            dx = 0.d0
            do j = 1, nno
                dx = dx+dff(i, j)*tabco(ndim*(j-1)+k)
            end do
            jac(2, i) = jac(2, i)+dx*nor(k)
        end do
    end do
!
! --- CALCUL DE R3,DR2 EN KSI EN FONTION DE LA FACE A APPARTENIR
! ---           R3 : LA COTE DE L'ELEMENT ENFANT
!
    if (ndime .eq. 3) then
! --- RECUPERATION DES COORDONNEES DANS L'ESPACE DE REFERENCE DES
! --- NOEUDS DE L'ELEMENT PARENT
        refcoo = 0.d0
        call elrfno(elp, nodeCoor=refcoo)
!
! ---  CALCUL DES VECTEURS N(1)N(2) ET N(1)N(3)
        v1(:) = 0.d0
        v2(:) = 0.d0
        do i = 1, ndim
            v1(i) = refcoo(i, n(2))-refcoo(i, n(1))
            v2(i) = refcoo(i, n(3))-refcoo(i, n(1))
        end do
!
! calcul de la NORMALE a la face
        call provec(v1, v2, nf)
!
        do i = 1, ndim
            r(3) = r(3)+nf(i)*(ksi(i)-refcoo(i, n(1)))
            jac(3, i) = nf(i)
        end do
    end if
!
    det = 0.d0
    inv(:, :) = 0.d0
! --- CALCUL L'INVERSE DE LA MATRICE JACOBIENNE
    call matinv('S', ndime, jac, inv, det)
!
    delta(:) = 0.d0
! --- CALCUL DES QUANTITES A ANNULER
!     CALCUL DE DELTAS
    do i = 1, ndime
        do j = 1, ndime
            delta(i) = delta(i)+inv(i, j)*r(j)
        end do
    end do
!
end subroutine
