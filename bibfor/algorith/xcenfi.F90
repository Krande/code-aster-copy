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
subroutine xcenfi(elrefp, ndim, ndime, nno, geom, &
                  lsn, pinref, pmiref, cenref, cenfi, &
                  nn, exit, jonc, num)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/reerel.h"
#include "asterfort/xcedge.h"
#include "asterfort/xelrex.h"
#include "asterfort/xnewto.h"
#include "asterfort/xnormv.h"
    integer(kind=8) :: ndim, ndime, nno, nn(4), exit(2)
    integer(kind=8), intent(in), optional :: num(8)
    character(len=8) :: elrefp
    real(kind=8) :: lsn(*), geom(*), pinref(*), pmiref(*)
    real(kind=8) :: cenfi(ndim), cenref(ndime)
    aster_logical :: jonc
!                      TROUVER LES COORDONNES DU PT MILIEU
!                      D'UNE FACE QUADRANGLE
!
!     ENTREE
!       ELP     : TYPE DE L'ELEMENT
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       PTINT  : COORDONNÉES DES POINTS D'INTERSECTION
!       JTABCO  : ADRESSE DES COORDONNEES DES NOEUDS DE L'ELEMENT
!       JTABLS  : ADRESSE DES LSN DES NOEUDS DE L'ELEMENT
!       NNO     : NOMBRE DE NOEUX DE L'ELEMENT
!       NN      : NUMERO DES NOEUDS DU SOUS-TETRA
!       NUM     : NUMERO DES NOEUDS DE LA FACE QUADRANGLE
!       JONC    : L'ELEMENT CONTIENT-IL UNE JONCTION
!     SORTIE
!       CENFI   : COORDONNES DU PT MILIEU AU CENTRE DE LA FISSURE
!     ----------------------------------------------------------------
!
    real(kind=8) :: epsmax, rbid, crit, maxi, x(81), dekker(4*ndime)
    real(kind=8) :: dff(3, 27), gradls(ndime)
    real(kind=8) :: ptxx(2*ndime), ksi(ndime), tole, xmi(ndime)
    integer(kind=8) :: ibid, itemax, i, n(3), j
    integer(kind=8) :: pi1, pi2, pi3, pi4, m12, m13, m24, m34
    character(len=6) :: name
    character(len=3) :: edge
    aster_logical :: courbe
    parameter(tole=1.d-2)
!
! --------------------------------------------------------------------
!
!
    itemax = 100
    epsmax = 1.d-9
    name = 'XCENFI'
    if (present(num)) then
        pi1 = num(1)
        pi2 = num(2)
        pi3 = num(3)
        pi4 = num(4)
        m12 = num(5)
        m13 = num(6)
        m34 = num(7)
        m24 = num(8)
    else
!  CONFERE XSTUDO
        pi1 = 1
        pi2 = 2
        pi3 = 3
        pi4 = 4
        m12 = 9
        m13 = 12
        m34 = 11
        m24 = 10
    end if
!
    ASSERT(ndime .eq. 3)
!
!   CALCUL D UN POINT DE DEPART POUR LE NEWTON
!===============================================================
!   ON DEFINI UN CRITERE SIMPLE POUR PRENDRE EN COMPTE LA COURBURE
!   CE CRITERE EST EFFICACE QUAND L INTERSECTION DE SURFACE DE L ISO ZERO
!   ET DU SOUS TETRAEDRE FORME UN <<POLYEDRE>> NON CONVEXE
!   SITUATION PLUS PROBABLE LORS DE LA DECOUPE D UN TETRATRAEDRE
!   EN REVANCHE,
!   EN RAFINANT LE MAILLAGE, LA SURFACE DE LA LEVEL-SET DEVIENT PLANE DANS
!   LES SOUS-TETRAS, LE CRITERE CI DESSOUS N A PLUS D INCIDENCE SUR LE CALCUL
!
    courbe = .false.
    maxi = 0.d0
    edge = ""
!   ARETE I1-I2
    call xcedge(ndime, pinref, pi1, pi2, pmiref, &
                m12, crit)
    if (crit .gt. maxi) then
        maxi = crit
        edge = "A12"
    end if
!   ARETE I2-I4
    call xcedge(ndime, pinref, pi2, pi4, pmiref, &
                m24, crit)
    if (crit .gt. maxi) then
        maxi = crit
        edge = "A24"
    end if
!   ARETE I3-I4
    call xcedge(ndime, pinref, pi3, pi4, pmiref, &
                m34, crit)
    if (crit .gt. maxi) then
        maxi = crit
        edge = "A34"
    end if
!   ARETE I1-I3
    call xcedge(ndime, pinref, pi1, pi3, pmiref, &
                m13, crit)
    if (crit .gt. maxi) then
        maxi = crit
        edge = "A13"
    end if
!
    if (maxi .gt. tole) courbe = .true.
!
    if (.not. courbe) then
        do i = 1, ndime
            ptxx(i+ndime) = (pinref(ndime*(pi1-1)+i)+pinref(ndime*(pi4-1)+i))/2.d0
        end do
    else
        do i = 1, ndime
            if (edge .eq. "A12" .or. edge .eq. "A34") then
                ptxx(i+ndime) = (pmiref(ndime*(m12-1)+i)+pmiref(ndime*(m34-1)+i))/2.d0
            else if (edge .eq. "A13" .or. edge .eq. "A24") then
                ptxx(i+ndime) = (pmiref(ndime*(m13-1)+i)+pmiref(ndime*(m24-1)+i))/2.d0
            else
                ASSERT(.false.)
            end if
        end do
    end if
!
!    CALCUL DE LA DIRECTION DE RECHERCHE
!    ON CHOISIT LE GRADIENT DE LA LSN AU POINT DE DEPART
!
    do j = 1, ndime
        xmi(j) = ptxx(ndime+j)
        ptxx(j) = 0.d0
    end do
    call elrfdf(elrefp, xmi, dff, nno, ndim)
!
    gradls(:) = 0.d0
    do i = 1, nno
        do j = 1, ndime
            gradls(j) = gradls(j)+dff(j, i)*lsn(i)
        end do
    end do
    call xnormv(ndime, gradls, rbid)
    do j = 1, ndime
        ptxx(j) = gradls(j)
    end do
!
!    ON RENSEIGNE LES NOEUDS DU SOUS TETRA POUR LA METHODE DE DEKKER
    dekker(:) = 0.d0
    call xelrex(elrefp, nno, x)
    do j = 1, ndime
        dekker(j) = x(ndime*(nn(1)-1)+j)
        dekker(j+ndime) = x(ndime*(nn(2)-1)+j)
        dekker(j+2*ndime) = x(ndime*(nn(3)-1)+j)
        dekker(j+3*ndime) = x(ndime*(nn(4)-1)+j)
    end do
!!!!!ATTENTION INITIALISATION DU NEWTON:
    ksi(:) = 0.d0
    if (jonc) then
        call xnewto(elrefp, name, n, ndime, ptxx, &
                    ndim, geom, lsn, ibid, ibid, &
                    itemax, epsmax, ksi)
    else
        call xnewto(elrefp, name, n, ndime, ptxx, &
                    ndim, geom, lsn, ibid, ibid, &
                    itemax, epsmax, ksi, exit, dekker)
    end if
!
    do i = 1, ndime
        cenref(i) = ksi(1)*ptxx(i)+ptxx(i+ndime)
    end do
!
! --- COORDONNES DU POINT DANS L'ELEMENT REEL
    call reerel(elrefp, nno, ndim, geom, cenref, &
                cenfi)
!
end subroutine
