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
subroutine bptobg(m, n, p)
    implicit none
    real(kind=8) :: m(6), n(6), p(3, 3)
! ----------------------------------------------------------------------
! CHANGEMENT DE BASE DE LA BASE PROPRE A LA BASE INITIALE DONT LES
!    VECTEURS PROPRES SONT RANGES DANS P POUR DES MATRICES SYMETRIQUES
!
! IN  M       : MATRICE DANS BP
! IN  P       : MATRICE DE PASSAGE B->BP
! OUT N       : MATRICE M DANS B
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, t(3, 3)
    real(kind=8) :: temp
    n(1) = 0.d0
    n(2) = 0.d0
    n(3) = 0.d0
    n(4) = 0.d0
    n(5) = 0.d0
    n(6) = 0.d0
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
    do i = 1, 3
        do j = 1, 3
            temp = m(t(i, j))
            n(1) = n(1)+p(1, i)*p(1, j)*temp
            n(2) = n(2)+p(2, i)*p(2, j)*temp
            n(3) = n(3)+p(3, i)*p(3, j)*temp
            n(4) = n(4)+p(1, i)*p(2, j)*temp
            n(5) = n(5)+p(1, i)*p(3, j)*temp
            n(6) = n(6)+p(2, i)*p(3, j)*temp
        end do
    end do
end subroutine
