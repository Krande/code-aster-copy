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
subroutine cmatve(mat, vectin, vectou, ndim)
    implicit none
!
!***********************************************************************
!    B. GUIGON    P. RICHARD                   DATE 06/04/92
!-----------------------------------------------------------------------
!  BUT:  < PRODUIT MATRICE VECTEUR COMPLEXE >
!
!   CETTE ROUTINE CALCULE LE PRODUIT D'UNE MATRICE COMPLEXE PAR UNE
! VECTEUR COMPLEXE
!
!-----------------------------------------------------------------------
!
! MAT      /I/: MATRICE COMPLEXE
! VECTIN   /I/: VECTEUR COMPLEXE D'ENTREE
! VECTOUT  /I/: VECTEUR COMPLXE DE SORTIE
! NDIM     /I/: DIMENSION DU VECTEUR ET DE LA MATRICE
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ndim
    integer(kind=8) :: i, j
    complex(kind=8) :: mat(ndim, ndim), vectin(ndim), vectou(ndim)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    do i = 1, ndim
        vectou(i) = dcmplx(0.d0, 0.d0)
        do j = 1, ndim
            vectou(i) = vectou(i)+mat(i, j)*vectin(j)
        end do
    end do
    goto 999
!
999 continue
end subroutine
