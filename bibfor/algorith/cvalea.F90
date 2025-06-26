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
subroutine cvalea(ndim, cmod, ndimax, nbmod)
    implicit none
!
!***********************************************************************
!    B. GUIGON   P. RICHARD                    DATE 06/04/92
!-----------------------------------------------------------------------
!  BUT:  < COMPLEXE VECTEUR ALEATOIRE >
!
!   CETTE ROUTINE INITIALISE UN ENSEMBLE DE VECTEURS STOCKES EN COLONNE
!   DANS UNE MATRICE COMPLEXE A DES VALEUR ALEATOIRE COMPLEXE DE MODULE
!   UNITAIRE
!
!-----------------------------------------------------------------------
!
! NDIM     /I/: DIMENSION EFFICACE DES VECTEURS
! CMOD     /O/: MODES PROPRES COMPLEXES SOLUTIONS
! NDIM     /I/: DIMENSION MAX DES VECTEURS DES VECTEURS
! NBMOD    /I/: NOMBRE DE VECTEURS
!
!-----------------------------------------------------------------------
!
#include "asterfort/ggubs.h"
    integer(kind=8) :: ndimax
    integer(kind=8) :: ndim, nbmod
    complex(kind=8) :: cmod(ndimax, nbmod)
    real(kind=8) :: r(2), dseed
    integer(kind=8) :: i, k
!
!-----------------------------------------------------------------------
    dseed = 24331.d0
    do k = 1, nbmod
        do i = 1, ndim
            call ggubs(dseed, 2, r)
            cmod(i, k) = dcmplx(r(1), r(2))
        end do
    end do
!
end subroutine
