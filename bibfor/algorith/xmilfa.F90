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
subroutine xmilfa(elrefp, ndim, ndime, geom, cnset, &
                  nnose, it, ainter, ip1, ip2, &
                  pm2, typma, pinref, pmiref, ksi, &
                  milfa, pintt, pmitt)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/conare.h"
#include "asterfort/reeref.h"
#include "asterfort/reerel.h"
#include "asterfort/xelrex.h"
#include "asterfort/xnormv.h"
#include "asterfort/xxmmvd.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: ip1, ip2, pm2, cnset(*), nnose, it, ndim, ndime
    real(kind=8) :: pinref(*), geom(*), milfa(ndim), ainter(*)
    real(kind=8) :: pmiref(*), ksi(ndime), pintt(*), pmitt(*)
    character(len=8) :: elrefp, typma
!                      TROUVER LES COORDONNES DU PT MILIEU ENTRE LES
!                      DEUX POINTS D'INTERSECTION
!
!     ENTREE
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       N       : NUMERO PARENT DES TROIS NOEUDS DE FACE
!       TYPMA   : TYPE DE LA MAILLE (TYPE_MAILLE)
!       IP1     : PREMIER IP DANS DANS UN FACE
!       IP2     : DEUXIEME IP DANS UN FACE
!       PM1A    : PREMIER PM1 DANS UN FACE
!       PM1B    : DEUXIEME PM1 DANS UN FACE
!       PM2     : PM2 DANS UN FACE
!       PINTER  : TABLEAU DES COORDONNÉES DES POINTS D'INTERSECTION
!       TABCO   : TABLEAU DES COORDONNÉES DES POINTS
!       PMILIE  : TABLEAU DES COORDONNÉES DES POINTS MILIEUX EXISTANT
!       AINTER  : INFOS ARETE ASSOCIEE AU POINT D'INTERSECTION
!       PINTT   : COORDONNEES REELES DES POINTS D'INTERSECTION
!       PMITT   : COORDONNEES REELES DES POINTS MILIEUX
!     SORTIE
!       MILFA   : COORDONNES DU TROISIME TYPE DE PM
!     ----------------------------------------------------------------
!
    integer(kind=8) :: a1, a2, a, b, d, ib, ar(12, 3), nbar, ia, id
    integer(kind=8) :: i, j, zxain, nno
    real(kind=8) :: xref(81), ptb(ndime), ptd(ndime), newpt(ndim)
    real(kind=8) :: pta(ndime), cosu, cosv, cosw
    real(kind=8) :: ff(27), t1(ndime), t2(ndime), sinu, rbid, t3(ndime)
    aster_logical :: courbe
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------
    zxain = xxmmvd('ZXAIN')
!   IDENTIFICATION DES NOEUDS DE LA FACE QUADRANGLE DANS LE SOUS TETRA
    call conare(typma, ar, nbar)
    a1 = nint(ainter(zxain*(ip1-1)+1))
    a2 = nint(ainter(zxain*(ip2-1)+1))
!
    a = 0
    b = 0
    d = 0
!
    do i = 1, 2
        do j = 1, 2
            if (ar(a1, i) .eq. ar(a2, j)) then
                a = ar(a1, 3-i)
                b = ar(a2, 3-j)
            end if
        end do
    end do
    do i = 1, nbar
        do j = 1, 2
            if ((ar(i, j) .eq. a) .and. (ar(i, 3-j) .eq. b)) d = ar(i, 3)
        end do
    end do
    ASSERT((a*b*d) .gt. 0)
!   INDICE CORRECPONDANT DANS L ELEMENT PARENT
    ia = cnset(nnose*(it-1)+a)
    ib = cnset(nnose*(it-1)+b)
    id = cnset(nnose*(it-1)+d)
!
    call xelrex(elrefp, nno, xref)
!
    if (ib .lt. 1000) then
        do j = 1, ndime
            ptb(j) = xref(ndime*(ib-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pintt(ndim*(ib-1001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    ptb, ff)
    end if
!
    if (id .lt. 2000) then
        do j = 1, ndime
            ptd(j) = xref(ndime*(id-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pmitt(ndim*(id-2001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    ptd, ff)
    end if
!
    if (ia .lt. 1000) then
        do j = 1, ndime
            pta(j) = xref(ndime*(ia-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pintt(ndim*(ia-1001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    pta, ff)
    end if
!
    do i = 1, ndime
        ksi(i) = (pinref(ndime*(ip1-1)+i)+ptb(i))/2.d0
    end do
! --- TEST SI LSN COURBE :
    courbe = .false.
    do i = 1, ndime
        t1(i) = ksi(i)-pinref(ndime*(ip1-1)+i)
        t2(i) = -1.5d0*pinref( &
                ndime*(ip1-1)+i)-5.d-1*pinref(ndime*(ip2-1)+i)+2.d0*pmiref(ndime*(pm2-1)+i)
        t3(i) = pta(i)-pinref(ndime*(ip1-1)+i)
    end do
    call xnormv(ndime, t1, rbid)
    call xnormv(ndime, t2, rbid)
    call xnormv(ndime, t3, rbid)
    b_n = to_blas_int(ndime)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    cosu = ddot(b_n, t1, b_incx, t2, b_incy)
    sinu = sqrt(1-cosu**2)
!   ON CHOISIT UNE CONVENTION DE SIGNE
    b_n = to_blas_int(ndime)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    cosv = ddot(b_n, t3, b_incx, t2, b_incy)
    b_n = to_blas_int(ndime)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    cosw = ddot(b_n, t3, b_incx, t1, b_incy)
    if (cosv .gt. cosw) sinu = -sinu
!
!   ON RAJOUTE UNE TOLE POUR EVITER DES DECOUPES TROP POURRIES
    if (sinu .lt. 1.d-3) courbe = .true.
!
    if (courbe) then
!   EN DEUXIEME APPROXIMATION: ON CHOISIT LE MILIEU DES "MILIEUX" PM2 ET D
!
        do i = 1, ndime
            ksi(i) = (pmiref(ndime*(pm2-1)+i)+ptd(i))/2.d0
        end do
    end if
!
! --- COORDONNES DU POINT DANS L'ELEMENT REEL
!
    do i = 1, ndime
        ASSERT(abs(ksi(i)) .le. 1.d0+1.d-12)
        if (ksi(i) .gt. 1.d0) ksi(i) = 1.d0
        if (ksi(i) .lt. -1.d0) ksi(i) = -1.d0
    end do
    call reerel(elrefp, nno, ndim, geom, ksi, &
                milfa)
end subroutine
