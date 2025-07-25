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
subroutine trndgl(nbx, vectn, vectpt, deplg, depll, &
                  rotfic)
    implicit none
!
#include "jeveux.h"
    integer(kind=8) :: nbx
    real(kind=8) :: vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: deplg(*), depll(*), t(3, 3), rotfic(*)
    integer(kind=8) :: i, i1, i2, ib, j
!-----------------------------------------------------------------------
!
    do ib = 1, nbx
!
!     RESTITUTION DE LA MATRICE DE PASSAGE
!
        do i = 1, 2
            do j = 1, 3
                t(i, j) = vectpt(ib, i, j)
            end do
        end do
        t(3, 1) = vectn(ib, 1)
        t(3, 2) = vectn(ib, 2)
        t(3, 3) = vectn(ib, 3)
!
        i1 = 5*(ib-1)
        i2 = 6*(ib-1)
!
!     LES TERMES DE TRANSLATION
!
        if (ib .le. nbx-1) then
            depll(i1+1) = deplg(i2+1)
            depll(i1+2) = deplg(i2+2)
            depll(i1+3) = deplg(i2+3)
!
!     LES TERMES DE ROTATION (2 SEULEMENT)
!
            depll(i1+4) = t(1, 1)*deplg(i2+4)+t(1, 2)*deplg(i2+5)+t(1, 3)* &
                          deplg(i2+6)
!
            depll(i1+5) = t(2, 1)*deplg(i2+4)+t(2, 2)*deplg(i2+5)+t(2, 3)* &
                          deplg(i2+6)
!
            rotfic(ib) = t(3, 1)*deplg(i2+4)+t(3, 2)*deplg(i2+5)+t(3, 3)* &
                         deplg(i2+6)
        else
            depll(i1+1) = t(1, 1)*deplg(i2+1)+t(1, 2)*deplg(i2+2)+t(1, 3)* &
                          deplg(i2+3)
!
            depll(i1+2) = t(2, 1)*deplg(i2+1)+t(2, 2)*deplg(i2+2)+t(2, 3)* &
                          deplg(i2+3)
!
            rotfic(ib) = t(3, 1)*deplg(i2+1)+t(3, 2)*deplg(i2+2)+t(3, 3)* &
                         deplg(i2+3)
        end if
!
    end do
!
end subroutine
