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

subroutine vecini(n, s, x)
    implicit none
!       INITIALISATION DU VECTEUR   X = S
!       IN  S      :  SCALAIRE
!       IN  N      :  DIMENSION DE X
!       OUT X      :  VECTEUR RESULTAT
!       POUR TOUS LES TYPES DE DONNEES VOIR AUSSI VECINI, VECINT, VECINK
!       ET VECINC.
!       ----------------------------------------------------------------
    integer(kind=8) :: n, i
    real(kind=8) :: x(n), s
    do i = 1, n
        x(i) = s
    end do
end subroutine
