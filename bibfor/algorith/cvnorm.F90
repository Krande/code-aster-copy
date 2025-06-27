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
subroutine cvnorm(mat, vect, ndim, iretou)
    implicit none
!
!***********************************************************************
!    B. GUIGON    P. RICHARD                   DATE 06/04/92
!-----------------------------------------------------------------------
!  BUT:  < NORME VECTEUR >
!
!   CETTE ROUTINE NORME UN VECTEUR COMPLEXE PAR RAPPORT A UNE MATRICE
!   COMPLEXE
!
!-----------------------------------------------------------------------
!
! MAT      /I/: MATRICE COMPLEXE DEFINISSANT LE PRODUIT SCALAIRE
! VECT     /M/: VECTEUR COMPLEXE A NORMER
! NDIM     /I/: DIMENSION DU VECTEUR ET DE LA MATRICE
! IRETOU
!-----------------------------------------------------------------------
!
#include "asterfort/sesqui.h"
    integer(kind=8) :: ndim
    complex(kind=8) :: mat(*), vect(ndim)
    real(kind=8) :: zero
    integer(kind=8) :: i, iretou
    complex(kind=8) :: normec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    data zero/0.d0/
!
!-----------------------------------------------------------------------
!
    call sesqui(mat, vect, ndim, normec)
    iretou = 0
    if (abs(normec) .eq. zero) then
        iretou = 1
        goto 999
    end if
    normec = dcmplx(sqrt(abs(dble(normec))), 0.d0)
    do i = 1, ndim
        vect(i) = vect(i)/normec
    end do
999 continue
end subroutine
