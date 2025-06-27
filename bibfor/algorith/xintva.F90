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
subroutine xintva(name, dekker, ptxx, ndime, intinf, &
                  intsup)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/provec.h"
#include "asterfort/xnormv.h"
#include "blas/ddot.h"
    integer(kind=8) :: ndime
    character(len=6) :: name
    real(kind=8) :: ptxx(*), intinf, intsup, dekker(4*ndime)
!
!         RECHERCHE DE L'INTERVALE DE RECHERCHE POUR LA RESOLUTION
!                   DU PROBLEME SCALAIRE DANS XMIFIS
!
!     ENTREE
!       N       : INDICE DES NOEUDS SOMMETS DE LA FACE TRIA DANS
!                 LAQUELLE ON EFFECTUE LA RECHERCHE
!       PTXX    : PONT DE DEPART ET DIRECTION DE RECHERCHE
!     SORTIE
!       INTINF  : BORNE INFERIEURE DE L'INTERVALLE DE RECHERCHE
!       INTSUP  : BORNE SUPERIEURE DE L'INTERVALLE DE RECHERCHE
!     --------------------------------------------------------------
!
    real(kind=8) :: a(ndime), b(ndime), c(ndime), ab(ndime), bc(ndime)
    real(kind=8) :: ca(ndime), ptini(ndime), k(ndime), det, alpha, kappa
    real(kind=8) :: ka(ndime), kb(ndime), kc(ndime), b1, c1, c2, ad(ndime)
    real(kind=8) :: norm_ab, norm_bc, norm_ca, eps, tole, d(ndime), bd(ndime)
    real(kind=8) :: norm_n, n(ndime), norm_ad, norm_bd
    integer(kind=8) :: j, cpt
    blas_int :: b_incx, b_incy, b_n
    parameter(eps=1.d-12)
    parameter(tole=5.d-7)
!
! ------------------------------------------------------------------
!
    call jemarq()
!
    if (name .eq. 'XINTER') then
        intinf = 0.d0
        intsup = 1.d0
        goto 99
    end if
!
    intinf = 0.d0
    intsup = 0.d0
    cpt = 0
!   RECUPERATION DES COORDONNEES DE REFERENCE DES NOEUDS SOMMETS DU SOUS ELEMENT
    do j = 1, ndime
        a(j) = dekker(j)
        b(j) = dekker(ndime+j)
        c(j) = dekker(2*ndime+j)
        d(j) = dekker(3*ndime+j)
    end do
!
!   RECUPERATION DES ARETES DU SOUS ELEMENT
    do j = 1, ndime
        ab(j) = b(j)-a(j)
        bc(j) = c(j)-b(j)
        ca(j) = a(j)-c(j)
        ad(j) = d(j)-a(j)
        bd(j) = d(j)-b(j)
    end do
!
    call xnormv(ndime, ab, norm_ab)
    call xnormv(ndime, bc, norm_bc)
    call xnormv(ndime, ca, norm_ca)
    call xnormv(ndime, ad, norm_ad)
    call xnormv(ndime, bd, norm_bd)
!
!   RECUPERACTION DU POINT DE DEPART ET DE LA DIRECTION DE RECHERCHE
    do j = 1, ndime
        ptini(j) = ptxx(ndime+j)
        k(j) = ptxx(j)
    end do
!
    if (name .eq. 'XMIFIS') then
!   ON CHERCHE LES INTERSECTIONS DE LA DIRECTION DE RECHERCHE AVEC LES ARETES DE LA FACE TRIA
!
!   INTERSECTION SUR L'ARETE AB (ON PROJETTE LE SYSTEME D'EQUATION DANS LE REPERE DE LA FACE)
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        b1 = ddot(b_n, k, b_incx, ab, b_incy)
        do j = 1, ndime
            ka(j) = a(j)-ptini(j)
        end do
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        c1 = ddot(b_n, k, b_incx, ka, b_incy)
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        c2 = ddot(b_n, ka, b_incx, ab, b_incy)
        det = -1.d0+b1**2
        if (abs(det) .gt. eps) then
            kappa = (-c1+c2*b1)/det
            alpha = (c2-b1*c1)/det
            alpha = alpha/norm_ab
            if (alpha .le. (1.d0+tole) .and. alpha .ge. -tole) then
                if (kappa .le. 0.d0) intinf = kappa
                if (kappa .ge. 0.d0) intsup = kappa
                cpt = cpt+1
            end if
        end if
!
!   INTERSECTION SUR L'ARETE BC
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        b1 = ddot(b_n, k, b_incx, bc, b_incy)
        do j = 1, ndime
            kb(j) = b(j)-ptini(j)
        end do
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        c1 = ddot(b_n, k, b_incx, kb, b_incy)
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        c2 = ddot(b_n, kb, b_incx, bc, b_incy)
        det = -1.d0+b1**2
        if (abs(det) .gt. eps) then
            kappa = (-c1+c2*b1)/det
            alpha = (c2-b1*c1)/det
            alpha = alpha/norm_bc
            if (alpha .le. (1.d0+tole) .and. alpha .ge. -tole) then
                if (kappa .le. 0.d0) intinf = kappa
                if (kappa .ge. 0.d0) intsup = kappa
                cpt = cpt+1
            end if
        end if
!
!   INTERSECTION SUR L'ARETE CA
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        b1 = ddot(b_n, k, b_incx, ca, b_incy)
        do j = 1, ndime
            kc(j) = c(j)-ptini(j)
        end do
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        c1 = ddot(b_n, k, b_incx, kc, b_incy)
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        c2 = ddot(b_n, kc, b_incx, ca, b_incy)
        det = -1.d0+b1**2
        if (abs(det) .gt. eps) then
            kappa = (-c1+c2*b1)/det
            alpha = (c2-b1*c1)/det
            alpha = alpha/norm_ca
            if (alpha .le. (1.d0+tole) .and. alpha .ge. -tole) then
                if (kappa .le. 0.d0) intinf = kappa
                if (kappa .ge. 0.d0) intsup = kappa
                cpt = cpt+1
            end if
        end if
    else if (name .eq. 'XCENFI') then
!   ON CHERCHE LES INTERSECTIONS DE LA DIRECTION DE RECHERCHE AVEC LES FACES DU TETRA
        intinf = -10.d0
        intsup = 10.d0
!   INTERSECTION SUR LA FACE ABC
!      RECHERCHE D'UN VECTEUR UNITAIRE NORMAL A LA FACE
        call provec(ab, bc, n)
        call xnormv(ndime, n, norm_n)
!      SI LA DIRECTION DE RECHERCHE EST PARALLELE A LA FACE ON SORT
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        if (abs(ddot(b_n, n, b_incx, k, b_incy)) .ge. tole) then
            do j = 1, ndime
                ka(j) = a(j)-ptini(j)
            end do
            b_n = to_blas_int(ndime)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            alpha = ddot(b_n, ka, b_incx, n, b_incy)/ddot(b_n, n, b_incx, k, b_incy)
            if (alpha .gt. intinf .and. alpha .lt. 0.d0) then
                intinf = alpha
                cpt = cpt+1
            else if (alpha .lt. intsup .and. alpha .gt. 0.d0) then
                intsup = alpha
                cpt = cpt+1
            end if
        end if
!
!   INTERSECTION SUR LA FACE ADB
!      RECHERCHE D'UN VECTEUR UNITAIRE NORMAL A LA FACE
        call provec(ab, ad, n)
        call xnormv(ndime, n, norm_n)
!      SI LA DIRECTION DE RECHERCHE EST PARALLELE A LA FACE ON SORT
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        if (abs(ddot(b_n, n, b_incx, k, b_incy)) .ge. tole) then
            do j = 1, ndime
                ka(j) = a(j)-ptini(j)
            end do
            b_n = to_blas_int(ndime)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            alpha = ddot(b_n, ka, b_incx, n, b_incy)/ddot(b_n, n, b_incx, k, b_incy)
            if (alpha .gt. intinf .and. alpha .lt. 0.d0) then
                intinf = alpha
                cpt = cpt+1
            else if (alpha .lt. intsup .and. alpha .gt. 0.d0) then
                intsup = alpha
                cpt = cpt+1
            end if
        end if
!
!   INTERSECTION SUR LA FACE ACD
!      RECHERCHE D'UN VECTEUR UNITAIRE NORMAL A LA FACE
        call provec(ca, ad, n)
        call xnormv(ndime, n, norm_n)
!      SI LA DIRECTION DE RECHERCHE EST PARALLELE A LA FACE ON SORT
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        if (abs(ddot(b_n, n, b_incx, k, b_incy)) .ge. tole) then
            do j = 1, ndime
                ka(j) = a(j)-ptini(j)
            end do
            b_n = to_blas_int(ndime)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            alpha = ddot(b_n, ka, b_incx, n, b_incy)/ddot(b_n, n, b_incx, k, b_incy)
            if (alpha .gt. intinf .and. alpha .lt. 0.d0) then
                intinf = alpha
                cpt = cpt+1
            else if (alpha .lt. intsup .and. alpha .gt. 0.d0) then
                intsup = alpha
                cpt = cpt+1
            end if
        end if
!
!   INTERSECTION SUR LA FACE BCD
!      RECHERCHE D'UN VECTEUR UNITAIRE NORMAL A LA FACE
        call provec(bd, bc, n)
        call xnormv(ndime, n, norm_n)
!      SI LA DIRECTION DE RECHERCHE EST PARALLELE A LA FACE ON SORT
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        if (abs(ddot(b_n, n, b_incx, k, b_incy)) .ge. tole) then
            do j = 1, ndime
                kb(j) = b(j)-ptini(j)
            end do
            b_n = to_blas_int(ndime)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            alpha = ddot(b_n, kb, b_incx, n, b_incy)/ddot(b_n, n, b_incx, k, b_incy)
            if (alpha .gt. intinf .and. alpha .lt. 0.d0) then
                intinf = alpha
                cpt = cpt+1
            else if (alpha .lt. intsup .and. alpha .gt. 0.d0) then
                intsup = alpha
                cpt = cpt+1
            end if
        end if
    end if
!
    ASSERT(cpt .ge. 2)
!
99  continue
!
    call jedema()
end subroutine
