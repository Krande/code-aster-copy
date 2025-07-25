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
subroutine bcoude(igau, icou, isect, l, h, &
                  a, m, nno, ncou, nsect, &
                  ff, df1, df2, mmt, b)
    implicit none
#include "asterc/r8pi.h"
    real(kind=8) :: l, h, a, b(4, *), ff(*), df1(*), df2(*)
    integer(kind=8) :: nno, ncou, nsect, m, igau, icou, isect, mmt
!
! ......................................................................
!
!   BUT      : CALCUL DE LA MATRICE DE DEFORMATION B POUR LES ELEMENTS
!              TUYAU
!
!   ENTREES : <----
!           IGAU,ICOU,ISECT : LES INDICES DES POINTS D'INTEGRATION
!                RESPECTIVEMENT SUR LA LONGUEUR, DANS L'EPAISSEUR ET
!                SUR LA CIRCONFERENCE
!           NNO : NBRES DE NOEUDS
!           NCOU,NSECT : NBRES DE COUCHES, NBRES DE SECTEURS
!           FF : VALEURS DES FONCTIONS DE FORMES
!           DF1: VALEURS DES DERIVEES 1ERES DES FONCTIONS DE FORMES
!           DF2: VALEURS DES DERIVEES 2EMES DES FONCTIONS DE FORMES
!           M : NBRES DE MODES DE FOURIER
!           L : LONGEUR DE L'ELEMENT
!           H : SON EPAISSEUR
!           A : SON RAYON MOYEN
!   SORTIES : ---->
!           B : LA MATRICE DE DEFORMATIONS
! ......................................................................
!
    real(kind=8) :: pi, deuxpi, cosfi, sinfi, cosmfi, sinmfi
    real(kind=8) :: zeta, r, hk, dhk, d2hk, fi
    integer(kind=8) :: k, icolon, n, ibloc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    pi = r8pi()
    deuxpi = 2.d0*pi
!
    zeta = (icou-1)*h/(2.d0*ncou)-h/2.d0
    fi = (isect-1)*deuxpi/(2.d0*nsect)
    cosfi = cos(fi)
    sinfi = sin(fi)
    if (mmt .eq. 0) then
        r = a
    else
        r = a+zeta
    end if
!
!  REMLISSAGE DE LA MATRICE
!
    do k = 1, nno
!
        hk = ff(nno*(igau-1)+k)
        dhk = df1(nno*(igau-1)+k)*(2.d0/l)
        d2hk = df2(nno*(igau-1)+k)*(2.d0/l)*(2.d0/l)
!
! LE 2/L EST DU AU PASSAGE DE L'ELEMENT DE REFERENCE A L'ELEMENT
! REEL
!
        ibloc = (9+6*(m-1))*(k-1)
!
!  PARTIE POUTRE
!
!   1 ERE LIGNE
!
        b(1, ibloc+1) = dhk
        b(1, ibloc+2) = 0.d0
        b(1, ibloc+3) = 0.d0
        b(1, ibloc+4) = 0.d0
        b(1, ibloc+5) = -r*cosfi*dhk
        b(1, ibloc+6) = r*sinfi*dhk
!
!   2 EME LIGNE
!
        b(2, ibloc+1) = 0.d0
        b(2, ibloc+2) = 0.d0
        b(2, ibloc+3) = 0.d0
        b(2, ibloc+4) = 0.d0
        b(2, ibloc+5) = 0.d0
        b(2, ibloc+6) = 0.d0
!
!   3 EME LIGNE
!
        b(3, ibloc+1) = 0.d0
        b(3, ibloc+2) = -cosfi*dhk
        b(3, ibloc+3) = sinfi*dhk
        b(3, ibloc+4) = -r*dhk
        b(3, ibloc+5) = sinfi*hk
        b(3, ibloc+6) = cosfi*hk
!
!   4 EME LIGNE
!
        b(4, ibloc+1) = 0.d0
        b(4, ibloc+2) = -sinfi*dhk
        b(4, ibloc+3) = -cosfi*dhk
        b(4, ibloc+4) = 0.d0
        b(4, ibloc+5) = -cosfi*hk
        b(4, ibloc+6) = sinfi*hk
!
!  FIN PARTIE POUTRE ET DEBUT PARTIE SUPPLEMENTAIRE
!
!      PARTIE IN-PLANE
!
        do n = 2, m
            icolon = ibloc+6+3*(n-2)
            cosmfi = cos(n*fi)
            sinmfi = sin(n*fi)
!
            b(1, icolon+1) = dhk*cosmfi
            b(1, icolon+2) = 0.d0
            b(1, icolon+3) = -zeta*d2hk*cosmfi
!
            b(2, icolon+1) = 0.d0
            b(2, icolon+2) = (n/r)*hk*cosmfi*(1.d0+zeta/a)
            b(2, icolon+3) = (1.d0/r)*hk*cosmfi*(1.d0+zeta*n*n/a)
!
!
            b(3, icolon+1) = -(n/r)*hk*sinmfi
            b(3, icolon+2) = dhk*sinmfi*(1+zeta/a)
            b(3, icolon+3) = zeta*n*dhk*sinmfi*(1.d0/a+1.d0/r)
!
!
            b(4, icolon+1) = 0.d0
            b(4, icolon+2) = 0.d0
            b(4, icolon+3) = 0.d0
!
!
        end do
!
!  FIN PARTIE IN-PLANE DEBUT PARTIE OUT-OF-PLANE
!
!
        do n = 2, m
            icolon = ibloc+6+3*(m-1)+3*(n-2)
            cosmfi = cos(n*fi)
            sinmfi = sin(n*fi)
!
            b(1, icolon+1) = dhk*sinmfi
            b(1, icolon+2) = 0.d0
            b(1, icolon+3) = -zeta*d2hk*sinmfi
!
            b(2, icolon+1) = 0.d0
            b(2, icolon+2) = -(n/r)*hk*sinmfi*(1.d0+zeta/a)
            b(2, icolon+3) = (1.d0/r)*hk*sinmfi*(1.d0+zeta*n*n/a)
!
!
            b(3, icolon+1) = (n/r)*hk*cosmfi
            b(3, icolon+2) = dhk*cosmfi*(1.d0+zeta/a)
            b(3, icolon+3) = -zeta*n*dhk*cosmfi*(1.d0/a+1.d0/r)
!
!
            b(4, icolon+1) = 0.d0
            b(4, icolon+2) = 0.d0
            b(4, icolon+3) = 0.d0
!
!
        end do
!
!  FIN OUT-OF-PLANE DEBUT PARTIE GONFLEMENT
!
        icolon = ibloc+6*(m-1)+6
        b(1, icolon+1) = -zeta*d2hk
        b(2, icolon+1) = hk/r
        b(3, icolon+1) = 0.d0
        b(4, icolon+1) = 0.d0
!
!  FIN PARTIE GONFLEMENT
!
!  1ERS MODES EN W
!
        b(1, icolon+2) = -zeta*d2hk*cosfi
        b(2, icolon+2) = (2.d0/r)*hk*cosfi*(1.d0+zeta/a)
        b(3, icolon+2) = dhk*sinfi*(1.d0+2.d0*zeta/a+zeta/r)
        b(4, icolon+2) = 0.d0
!
        b(1, icolon+3) = -zeta*d2hk*sinfi
        b(2, icolon+3) = (2.d0/r)*hk*sinfi*(1.d0+zeta/a)
        b(3, icolon+3) = -dhk*cosfi*(1.d0+2.d0*zeta/a+zeta/r)
        b(4, icolon+3) = 0.d0
!
!
!  FIN PARTIE GONFLEMENT
!
! FIN DES 1ERS MODES EN W
!
    end do
!
!
! FIN REMPLISSAGE DE LA MATRICE
!
!
end subroutine
