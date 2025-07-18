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

subroutine dhrc_recup_mate(imate, compor, a0, c0, aa_t, ga_t, ab, gb, &
                           ac, gc, aa_c, ga_c, cstseu)
! person_in_charge: sebastien.fayolle at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/utmess.h"
#include "asterfort/rcadlv.h"
#include "asterfort/assert.h"
!
    character(len=16), intent(in) :: compor
    integer(kind=8), intent(in) :: imate
    real(kind=8), intent(out) :: a0(6, 6), c0(2, 2, 2)
    real(kind=8), intent(out) :: aa_t(6, 6, 2), ab(6, 2, 2), ac(2, 2, 2)
    real(kind=8), intent(out) :: ga_t(6, 6, 2), gb(6, 2, 2), gc(2, 2, 2)
    real(kind=8), intent(out) :: aa_c(6, 6, 2)
    real(kind=8), intent(out) :: ga_c(6, 6, 2)
    real(kind=8), intent(out) :: cstseu(6)
! ----------------------------------------------------------------------
!
! BUT : LECTURE DES PARAMETRES MATERIAU POUR LE MODELE DHRC
!
! IN:
!       IMATE   : ADRESSE DU MATERIAU
!       COMPOR  : COMPORTMENT
!       EP      : EPAISSEUR DE LA PLAQUE
! OUT:
!       A0      : PARAMETRE D ELASTICITE
!       C0      : PARAMETRE D ELASTICITE
!       AA_T    : PARAMETRE ALPHA POUR LE TENSEUR A EN TRACTION
!       AA_C    : PARAMETRE ALPHA POUR LE TENSEUR A EN COMPRESSION
!       AB      : PARAMETRE ALPHA POUR LE TENSEUR B
!       AC      : PARAMETRE ALPHA POUR LE TENSEUR C
!       GA_T    : PARAMETRE GAMMA POUR LE TENSEUR A EN TRACTION
!       GA_C    : PARAMETRE GAMMA POUR LE TENSEUR A EN COMPRESSION
!       GB      : PARAMETRE GAMMA POUR LE TENSEUR B
!       GC      : PARAMETRE GAMMA POUR LE TENSEUR C
!
!       CSTSEU  : PARAMETRES DE SEUILS
!            (1): POUR L'ENDOMMAGEMENT
!            (2): POUR LE GLISSEMENT
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, l, jadr, n1, icodre1
! ----------------------------------------------------------------------
!
    if ((.not. (compor(1:4) .eq. 'DHRC'))) then
        call utmess('F', 'ELEMENTS4_65', sk=compor)
    end if
!
!     -----------------------------------------------------------------
!     MATRICE A0(6,6)
!     -----------------------------------------------------------------
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'ELAS_DHRC', 'A0', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 21)
!
    l = 0
    do i = 1, 6
        do j = i, 6
            l = l+1
            a0(j, i) = zr(jadr-1+l)
            a0(i, j) = a0(j, i)
        end do
    end do
    ASSERT(l .eq. 21)
!
!     -----------------------------------------------------------------
!     MATRICE AA_C(6,6,1)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'AA_C', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 42)
!
    l = 0
    do i = 1, 6
        do j = i, 6
            l = l+1
            aa_c(j, i, 1) = zr(jadr-1+l)
            aa_c(i, j, 1) = aa_c(j, i, 1)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE AA_C(6,6,2)
!     -----------------------------------------------------------------
!
    do i = 1, 6
        do j = i, 6
            l = l+1
            aa_c(j, i, 2) = zr(jadr-1+l)
            aa_c(i, j, 2) = aa_c(j, i, 2)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE AA_T(6,6,1)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'AA_T', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 42)
!
    l = 0
    do i = 1, 6
        do j = i, 6
            l = l+1
            aa_t(j, i, 1) = zr(jadr-1+l)
            aa_t(i, j, 1) = aa_t(j, i, 1)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE AA_T(6,6,2)
!     -----------------------------------------------------------------
!
    do i = 1, 6
        do j = i, 6
            l = l+1
            aa_t(j, i, 2) = zr(jadr-1+l)
            aa_t(i, j, 2) = aa_t(j, i, 2)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE GA_C(6,6,1)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'GA_C', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 42)
!
    l = 0
    do i = 1, 6
        do j = i, 6
            l = l+1
            ga_c(j, i, 1) = zr(jadr-1+l)
            ga_c(i, j, 1) = ga_c(j, i, 1)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE GA_C(6,6,2)
!     -----------------------------------------------------------------
!
    do i = 1, 6
        do j = i, 6
            l = l+1
            ga_c(j, i, 2) = zr(jadr-1+l)
            ga_c(i, j, 2) = ga_c(j, i, 2)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE GA_T(6,6,1)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'GA_T', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 42)
!
    l = 0
    do i = 1, 6
        do j = i, 6
            l = l+1
            ga_t(j, i, 1) = zr(jadr-1+l)
            ga_t(i, j, 1) = ga_t(j, i, 1)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE GA_T(6,6,2)
!     -----------------------------------------------------------------
!
    do i = 1, 6
        do j = i, 6
            l = l+1
            ga_t(j, i, 2) = zr(jadr-1+l)
            ga_t(i, j, 2) = ga_t(j, i, 2)
        end do
    end do
!
!     -----------------------------------------------------------------
!     MATRICE AB(6,2,2)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'AB', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 24)
    ab = reshape(source=zr(jadr:jadr+24), shape=(/6, 2, 2/))
!
!     -----------------------------------------------------------------
!     MATRICE GB(6,2,2)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'GB', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 24)
    gb = reshape(source=zr(jadr:jadr+24), shape=(/6, 2, 2/))
!
!     -----------------------------------------------------------------
!     MATRICE C0(2,2,2)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'C0', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 8)
    c0 = reshape(source=zr(jadr:jadr+8), shape=(/2, 2, 2/))
!
!     -----------------------------------------------------------------
!     MATRICE AC(2,2,2)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'AC', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 8)
    ac = reshape(source=zr(jadr:jadr+8), shape=(/2, 2, 2/))
!
!     -----------------------------------------------------------------
!     MATRICE GC(2,2,2)
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'GC', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 8)
    gc = reshape(source=zr(jadr:jadr+8), shape=(/2, 2, 2/))
!
!     -----------------------------------------------------------------
!     SEUILS NYD
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'NYD', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 2)
    cstseu(1:2) = zr(jadr:jadr+1)
!
!     -----------------------------------------------------------------
!     SEUILS SCRIT
!     -----------------------------------------------------------------
!
    call rcadlv(' ', 1, 1, '+', imate, ' ', 'DHRC', 'SCRIT', &
                0, [' '], [0.d0], jadr, n1, icodre1, 1)
    ASSERT(icodre1 .eq. 0 .and. n1 .eq. 4)
    cstseu(3:6) = zr(jadr:jadr+3)
!
end subroutine
