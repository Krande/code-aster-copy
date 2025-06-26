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
subroutine vpzech(d, z, low, high, mm, &
                  neq, iz)
    implicit none
    integer(kind=8) :: low, high, mm, neq, iz
    real(kind=8) :: d(1), z(iz, 1)
!     MISE A L'ECHELLE (NORMALISATION) DE LA COLONNE Z PAR LA BONNE
!     VALEUR DE "D"
!     ------------------------------------------------------------------
!     SERT A L'EQUILIBRAGE D'UNE MATRICE POUR LE CALCUL DE SES VALEURS
!     ET VECTEURS PROPRES.
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE 321 (ROUTINE BALBAK)
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, ii, jj
    real(kind=8) :: s
!     ------------------------------------------------------------------
    do i = low, high
        s = d(i)
        do j = 1, mm
            z(i, j) = z(i, j)*s
        end do
    end do
!
!     --- REPERMUTER LES LIGNES SI CA A ETE FAIT DANS VPZBAL ---
    do ii = low-1, 1, -1
        jj = nint(d(ii))
        if (ii .ne. jj) then
            do j = 1, mm
                s = z(ii, j)
                z(ii, j) = z(jj, j)
                z(jj, j) = s
            end do
        end if
    end do
!
    do ii = high+1, neq
        jj = nint(d(ii))
        if (ii .ne. jj) then
            do j = 1, mm
                s = z(ii, j)
                z(ii, j) = z(jj, j)
                z(jj, j) = s
            end do
        end if
    end do
end subroutine
