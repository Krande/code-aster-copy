! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine eifore(ndim, axi, nno1, nno2, npg, &
                  wref, vff1, vff2, dffr2, geom, &
                  ang, iu, im, sigref, depref, &
                  vect)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/eicine.h"
#include "asterfort/r8inir.h"
    aster_logical :: axi
    integer :: ndim, nno1, nno2, npg, iu(3, 18), im(3, 9)
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno2)
    real(kind=8) :: wref(npg)
    real(kind=8) :: vect(2*nno1*ndim+nno2*ndim)
    real(kind=8) :: dffr2(ndim-1, nno2, npg), ang(*), sigref, depref
! ----------------------------------------------------------------------
!
!     REFE_FORC_NODA POUR LES ELEMENTS D'INTERFACE
!
! ----------------------------------------------------------------------
! IN  NDIM    : DIMENSION DES ELEMENTS
! IN  AXI     : .TRUE. SI AXISYMETRIE
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE E)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE E)
! IN  DFFR2   : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE L)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  WREF    : POIDS DES POINTS DE GAUSS DE REFERENCE
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  IU      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! IN  IM      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE
! IN  SIGREF  : VALEUR DE REFERENCE POUR LES CONTRAINTES
! IN  DEPREF  : VALEUR DE REFERENCE POUR LES SAUTS DE DEPLACEMENT
! OUT VECT    : FORCES INTERIEURES DE REFERENCE
! ----------------------------------------------------------------------
    integer :: g, n, i, k, kk
    real(kind=8) :: wg, b(3, 3, 18), t1
! ----------------------------------------------------------------------
!
    call r8inir(nno1*2*ndim+nno2*ndim, 0.d0, vect, 1)
!
    do g = 1, npg
!
        call eicine(ndim, axi, nno1, nno2, vff1(1, g), &
                    vff2(1, g), wref(g), dffr2(1, 1, g), geom, ang, &
                    wg, b)
!
!      VECTEUR FINT:U
        do n = 1, 2*nno1
            do i = 1, ndim
                kk = iu(i, n)
                t1 = 0
                do k = 1, ndim
                    t1 = t1+abs(b(k, i, n))*sigref/ndim
                end do
                vect(kk) = vect(kk)+wg*t1
            end do
        end do
!
!      VECTEUR FINT:M
        do n = 1, nno2
            do i = 1, ndim
                kk = im(i, n)
                t1 = abs(vff2(n, g))*depref
                vect(kk) = vect(kk)+wg*t1
            end do
        end do
!
    end do
end subroutine
