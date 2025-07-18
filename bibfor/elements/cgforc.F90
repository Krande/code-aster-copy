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
subroutine cgforc(ndim, nno1, nno2, npg, wref, &
                  vff1, dffr1, geom, mat, pesa, &
                  iu, a, tang, vect)
!
    implicit none
!
#include "asterfort/cgcine.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
!
    integer(kind=8) :: ndim, nno1, nno2, npg, mat, iu(3, 3)
    real(kind=8) :: vff1(nno1, npg), geom(ndim, nno1), wref(npg)
    real(kind=8) :: vect(nno1*(ndim+1)+nno2)
    real(kind=8) :: dffr1(nno1, npg)
    real(kind=8) :: a, tang(3, 3), pesa(4)
! ----------------------------------------------------------------------
!
!   CHAR_MECA_PESA_R POUR ELEMENT CABLE/GAINE
!
! ----------------------------------------------------------------------
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE L)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  WREF    : POIDS DES POINTS DE GAUSS DE REFERENCE
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  DFFR1   : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE U)
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  TANG    : TANGENTE AUX NOEUDS
! IN  MAT     : MATERIAU CODE
! IN  IU      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT GA
! IN  A       : SECTION DE LA BARRE
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! ----------------------------------------------------------------------
    integer(kind=8) :: nddl, g, n, i, kk, codres(1)
    real(kind=8) :: wg, b(4, 3), t1
    real(kind=8) :: rho(1), courb, l(3)
! ----------------------------------------------------------------------
!
!
    nddl = nno1*(ndim+1)+nno2
    call r8inir(nddl, 0.d0, vect, 1)
!
!
! - CALCUL POUR CHAQUE POINT DE GAUSS
!
    do g = 1, npg
!
!      CALCUL DES ELEMENTS GEOM DE L'EF AU POINT DE GAUSS CONSIDERE
!
        call rcvalb('RIGI', g, 1, '+', mat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, 'RHO', rho, codres, 1)
!
        call cgcine(ndim, nno1, vff1(1, g), wref(g), dffr1(1, g), &
                    geom, tang, wg, l, b, &
                    courb)
!
!        VECTEUR FINT:U ET UC
        do n = 1, nno1
            do i = 1, ndim
                kk = iu(i, n)
                t1 = vff1(n, g)*a*rho(1)*pesa(1)*pesa(i+1)
                vect(kk) = vect(kk)+wg*t1
            end do
        end do
!
!
!
    end do
!
end subroutine
