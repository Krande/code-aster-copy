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
subroutine upletr(ndim, mple, mcol)
    implicit none
    integer(kind=8) :: ndim
    real(kind=8) :: mple(ndim, ndim), mcol(*)
!    CE SOUS PROGRAMME STOCKE SOUS FORME PLEINE UNE MATRICE
!    TRIANGULAIRE SUPERIEURE
!
!
!    -------------------------------------------------------------------
!
! IN TYPE ! NOM    ! TABLEAU !             SIGNIFICATION
! IN -------------------------------------------------------------------
! IN  I   ! NDIM   !     -   ! TAILLE DE LA MATRICE
! IN  R*8 !  MCOL  !    -    ! MATRICE UNICOLONNE TRIANGULAIRE
! SUPERIEURE DE TAILLE NDIM*(NDIM+1)/2
! IN
! IN (+) REMARQUES :
!
! OUT TYPE ! NOM   ! TABLEAU !             SIGNIFICATION
! OUT ------------------------------------------------------------------
! OUT R*8 ! MPLE   !NDIM*NDIM! MATRICE STOCKEE SOUS FORME PLEINE
!
!
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j
!
!
! ---------------------------------------------------------------------
! PASSAGE DE MATRICE COLONNE VERS MAT PLEINE (COMPLETEE PAR ANTI
! SYMETRIE)
    do i = 1, ndim
        do j = 1, i
            mple(j, i) = mcol(int(i*(i-1)/2)+j)
            mple(i, j) = -mple(j, i)
        end do
    end do
!
end subroutine
