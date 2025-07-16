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
subroutine nmgvdn(ndim, nno1, nno2, iu, ia)
!
!
!
    implicit none
!
    integer(kind=8) :: ndim, nno1, nno2, iu(ndim*nno1), ia(nno2)
! ---------------------------------------------------------------------
!
!     POSITION DES INDICES POUR LES DEGRES DE LIBERTE
!
! IN  NDIM    : DIMENSION DES ELEMENTS
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE A)
! ---------------------------------------------------------------------
    integer(kind=8) :: n, i, os
! ---------------------------------------------------------------------
!
!
!      ELEMENT P1 - CONTINU
!
    do n = 1, nno2
        do i = 1, ndim
            iu(nno1*(i-1)+n) = i+(n-1)*(ndim+1)
        end do
        ia(n) = 1+ndim+(n-1)*(ndim+1)
    end do
    os = (1+ndim)*nno2
    do n = 1, nno1-nno2
        do i = 1, ndim
            iu(nno1*(i-1)+n+nno2) = i+(n-1)*ndim+os
        end do
    end do
!
!
end subroutine
