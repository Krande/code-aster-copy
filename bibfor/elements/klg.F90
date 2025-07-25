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
subroutine klg(nno, nbrddl, pgl, k)
! aslint: disable=W1306
    implicit none
!
! PASSAGE DE LA MATRICE K DU REPERE LOCAL AU REPERE GLOBAL
! POUR LES DDL DE POUTRE, ON LAISSE LES DDL DE COQUE DANS
! LE REPERE LOCAL.
! ENTREE        : NNO = NBRE DE NOEUDS
!                 PGL = MATRICE DE PASSAGE
! ENTREE-SORTIE : K   = MATRICE DE RIGIDITE
!
    integer(kind=8) :: i, j, l, nno, nbrddl, m
!JMP      PARAMETER          (NBRDDL=63)
    real(kind=8) :: k(nbrddl, nbrddl), p(nbrddl, nbrddl), pgl(3, 3)
    real(kind=8) :: ktemp(nbrddl, nbrddl)
!
!  INITIALISATION A L'IDENTITE DE LA MATRICE DE PASSAGE P
!
    do i = 1, nbrddl
        do j = 1, nbrddl
            if (i .eq. j) then
                p(i, j) = 1.d0
            else
                p(i, j) = 0.d0
            end if
        end do
    end do
!
!  REMPLISSAGE DES DE BLOC DE LA MATRICE P CORRESPONDANT AUX DDL
!  DE POUTRE (UX, UY, UZ, TETAX, TETAY, ET TETAZ) PAR LA MATRICE
!  DE PASSAGE (3*3) PGL.
!
    do l = 1, nno
        m = (l-1)*nbrddl/nno
        do i = 1, 3
            do j = 1, 3
                p(m+i, m+j) = pgl(i, j)
                p(m+3+i, m+3+j) = pgl(i, j)
            end do
        end do
    end do
!
! INITIALISATION A ZERO DE LA MATRICE KTEMP
!
    do i = 1, nbrddl
        do j = 1, nbrddl
            ktemp(i, j) = 0.d0
        end do
    end do
!
! CALCUL DE KTEMP = PRODUIT (TRANSPOSEE P) * K
!
    do i = 1, nbrddl
        do j = 1, nbrddl
            do l = 1, nbrddl
                ktemp(i, j) = ktemp(i, j)+p(l, i)*k(l, j)
            end do
        end do
    end do
!
!  INITIALISATION A ZERO DE LA MATRICE K
!
    do i = 1, nbrddl
        do j = 1, nbrddl
            k(i, j) = 0.d0
        end do
    end do
!
! CALCUL DE K = PRODUIT KTEMP * P = (TRANSPOSEE P) * K * P
!
    do i = 1, nbrddl
        do j = 1, nbrddl
            do l = 1, nbrddl
                k(i, j) = k(i, j)+ktemp(i, l)*p(l, j)
            end do
        end do
    end do
!
!
end subroutine
