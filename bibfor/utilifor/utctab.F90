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
subroutine utctab(raz, na, mb, mc, a, &
                  b, c, xab, ctab)
    implicit none
#include "asterfort/r8inir.h"
    integer(kind=8) :: na, mb, mc
    character(len=*) :: raz
    real(kind=8) :: a(na, na), b(na, mb), c(na, mc), xab(na, mb)
    real(kind=8) :: ctab(mc, mb)
!     ------------------------------------------------------------------
!     PRODUIT CT . A . B - A CARREE - B ET C  RECTANGULAIRE
!     ------------------------------------------------------------------
!IN   K4  RAZ  'ZERO' : ON FAIT CTAB = 0    + CT*A.B
!              'CUMU' : ON FAIT CTAB = CTAB + CT*A.B
!IN   I   NA   ORDRE DE A
!IN   I   MB   NB DE COLONNES DE B
!IN   I   MC   NB DE COLONNES DE C
!IN   R   A    MATRICE A           (NA,NA)
!IN   R   B    MATRICE B           (NA,MB)
!IN   R   C    MATRICE C           (NA,MC)
!IN   R   XAB  ZONE DE TRAVAIL XAB (NA,MB)
!OUT  R   CTAB PRODUIT CT . A . B  (MC,MB)
!     ------------------------------------------------------------------
    character(len=4) :: raz2
! --DEB
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
!-----------------------------------------------------------------------
    raz2 = raz
!
    call r8inir(na*mb, 0.0d0, xab, 1)
    do i = 1, na
        do k = 1, na
            do j = 1, mb
                xab(i, j) = xab(i, j)+a(i, k)*b(k, j)
            end do
        end do
    end do
!
    if (raz2 .eq. 'ZERO') call r8inir(mc*mb, 0.0d0, ctab, 1)
!
    do i = 1, mc
        do k = 1, na
            do j = 1, mb
                ctab(i, j) = ctab(i, j)+c(k, i)*xab(k, j)
            end do
        end do
    end do
end subroutine
