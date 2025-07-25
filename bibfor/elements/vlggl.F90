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
subroutine vlggl(nno, nbrddl, pgl, v, code, &
                 p, vtemp)
    implicit none
#include "asterfort/utmess.h"
!
! PASSAGE D'UN VECTEUR V DU REPERE GLOBAL AU REPERE LOCAL
! OU INVERSEMENT. ON AGIT UNIQUEMENT SUR LES DDL DE POUTRE,
! LES DDL DE COQUE RESTENT INCHANGES.
!
    integer(kind=8) :: i, j, l, nno, nbrddl, m
!JMP      PARAMETER          (NBRDDL=63)
    real(kind=8) :: v(nbrddl), p(nbrddl, nbrddl), pgl(3, 3)
    real(kind=8) :: vtemp(nbrddl)
    character(len=2) :: code
!  ENTREE :NNO  = NBRE DE NOEUDS
!          PGL  = MATRICE DE PASSAGE
!          CODE = GL POUR UN PASSAGE GLOBAL -> LOCAL
!                 LG POUR UN PASSAGE LOCAL  -> GLOBAL
! ENTREE-SORTIE : V
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
! INITIALISATION A ZERO DU VECTEUR VTEMP
!
    do i = 1, nbrddl
        vtemp(i) = 0.d0
    end do
!
!  CAS D'UN PASSAGE LOCAL -> GLOBAL
!
    if (code .eq. 'LG') then
!
! CALCUL DE VTEMP = PRODUIT (TRANSPOSEE P) * V
!
        do i = 1, nbrddl
            do l = 1, nbrddl
                vtemp(i) = vtemp(i)+p(l, i)*v(l)
            end do
        end do
!
    else if (code .eq. 'GL') then
!
! CALCUL DE VTEMP = P * V
!
        do i = 1, nbrddl
            do l = 1, nbrddl
                vtemp(i) = vtemp(i)+p(i, l)*v(l)
            end do
        end do
!
    else
        call utmess('F', 'ELEMENTS4_58', sk=code)
    end if
!
! STOCKAGE DE VTEMP DANS V
!
    do i = 1, nbrddl
        v(i) = vtemp(i)
    end do
!
!
end subroutine
