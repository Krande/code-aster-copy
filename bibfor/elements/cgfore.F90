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
subroutine cgfore(ndim, nno1, nno2, npg, wref, &
                  vff1, vff2, dffr1, a, geom, &
                  tang, iu, iuc, im, forref, &
                  sigref, depref, vect)
!
    implicit none
!
#include "asterfort/cgcine.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: ndim, nno1, nno2, npg, iu(3, 3), iuc(3), im(3)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno1), wref(npg)
    real(kind=8) :: vect(2*nno1*ndim+nno2*ndim)
    real(kind=8) :: dffr1(nno1, npg), tang(*), sigref, depref, forref, a
! ----------------------------------------------------------------------
!
!     REFE_FORC_NODA POUR LES ELEMENTS CABLE GAINE
!
! ----------------------------------------------------------------------
! IN  NDIM    : DIMENSION DES ELEMENTS
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE E)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE E)
! IN  DFFR2   : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE L)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  WREF    : POIDS DES POINTS DE GAUSS DE REFERENCE
! IN  A       : SECTION CABLE
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  IU      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! IN  IM      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE
! IN  SIGREF  : VALEUR DE REFERENCE POUR LES CONTRAINTES
! IN  DEPREF  : VALEUR DE REFERENCE POUR LES SAUTS DE DEPLACEMENT
! OUT VECT    : FORCES INTERIEURES DE REFERENCE
! ----------------------------------------------------------------------
    integer(kind=8) :: g, n, i, kk
    real(kind=8) :: wg, b(4, 3), t1, l(3), courb, sig1r, sig2r, sig3r
! ----------------------------------------------------------------------
!
    call r8inir(nno1*(ndim+1)+nno2, 0.d0, vect, 1)
!
    do g = 1, npg
!
        call cgcine(ndim, nno1, vff1(1, g), wref(g), dffr1(1, g), &
                    geom, tang, wg, l, b, &
                    courb)
!        VECTEUR FINT:U
        sig1r = forref
        sig2r = sigref*sqrt(a)
        sig3r = depref
        do n = 1, nno1
!         POUR EVITER LES 0 SUR LES DDLS DX DY DZ DE GAINE
            t1 = 0.d0
            do i = 1, ndim
                t1 = t1+abs(b(i, n)*sig1r)
            end do
            do i = 1, ndim
                kk = iu(i, n)
                vect(kk) = vect(kk)+wg*t1/3.d0
            end do
            kk = iuc(n)
            t1 = abs(b(4, n)*sig1r)+abs(l(n)*sig2r)
            vect(kk) = vect(kk)+wg*t1
        end do
!
!        VECTEUR FINT:M
        do n = 1, nno2
            kk = im(n)
            t1 = abs(vff2(n, g)*sig3r)
            vect(kk) = vect(kk)+wg*t1
        end do
!
!
    end do
!
end subroutine
