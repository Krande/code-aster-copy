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
subroutine arlt1d(mlv, ndim, ndml2, mcpln2)
!
!
!
    implicit none
#include "jeveux.h"
!
    integer(kind=8) :: ndml2, ndim
    real(kind=8) :: mlv(78)
!
    real(kind=8) :: mcpln2(2*ndim*ndml2, 2*ndim*ndml2)
    integer(kind=8) :: iaux, jaux, kaux
!
! ----------------------------------------------------------------------
!
! CALCUL DES MATRICES DE COUPLAGE ARLEQUIN
! OPTION ARLQ_MATR : CALCUL DES INTEGRALES DE COUPLAGE 1D - 1D
! ----------------------------------------------------------------------
!
! --- CALCUL DES TERMES DE COUPLAGE - MATRICE STOCKAGE LINEAIRE
!
    do jaux = 1, 2*ndim*ndml2
        do iaux = 1, 2*ndim*ndml2
            mcpln2(iaux, jaux) = 0.d0
        end do
    end do
    kaux = 0
    do iaux = 1, 2*ndim*ndml2
        do jaux = 1, iaux
            kaux = kaux+1
            mcpln2(iaux, jaux) = mcpln2(iaux, jaux)+mlv(kaux)
            mcpln2(jaux, iaux) = mcpln2(iaux, jaux)
        end do
    end do
!
!
end subroutine
