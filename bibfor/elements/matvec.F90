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
subroutine matvec(nordre, amat, nombv, v1, v2, &
                  vecres)
    implicit none
#include "asterfort/utmess.h"
    real(kind=8) :: amat(*), v1(*), v2(*), vecres(*)
    real(kind=8) :: vsom(9)
! ---------------------------------------------
!     BUT:   POUR LES ELEMENTS DE CABLE, CALCUL DU PRODUIT DE LA MATRICE
!            AMAT D'ORDRE NORDRE PAR:
!            . SI NOMBV=1, LE VECTEUR V1;
!            . SI NOMBV=2, LA SOMME DES VECTEURS V1 ET V2.
!            LE VECTEUR RESULTAT EST MIS DANS VECRES.
!     IN: NORDRE
!         AMAT
!         NOMBV
!         V1
!         V2
!     OUT: VECRES
! ---------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, nombv, nordre
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    zero = 0.d0
    if (nombv .eq. 1) then
        do i = 1, nordre
            vsom(i) = v1(i)
        end do
    else if (nombv .eq. 2) then
        do i = 1, nordre
            vsom(i) = v1(i)+v2(i)
        end do
    else
        call utmess('F', 'ELEMENTS2_34')
    end if
    k = 0
    do i = 1, nordre
        vecres(i) = zero
        do j = 1, nordre
            k = k+1
            vecres(i) = vecres(i)+amat(k)*vsom(j)
        end do
    end do
end subroutine
