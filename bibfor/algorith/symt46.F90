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
subroutine symt46(ftot, ftos)
!
    implicit none
!
    real(kind=8), intent(in)  :: ftot(3, 3, 3, 3)
    real(kind=8), intent(out) :: ftos(6, 6)
! --------------------------------------------------------------------------------------------------
!     SYMETRISE UN TENSEUR SUR 2 PREMIERS ET 2 DERNIERS INDICES
!     NOTATION DE VOIGT AVEC RAC2 : XX,YY,ZZ,RAC2*XY,RAC2*XZ,RAC2*YZ
!
! IN  FTOT    : TENSEUR(3,3,3,3) F_IJKL AVEC SYMETRIES MINEURES
!               F_IJKL=F_JIKL et F_IJKL=F_IJLK
! OUT FTOS    : TENSEUR (6,6)
! - ------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, l, ijk(3, 3)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
!
! -- tableau de correspondance
!
    ijk(1, 1) = 1
    ijk(2, 2) = 2
    ijk(3, 3) = 3
    ijk(1, 2) = 4
    ijk(2, 1) = 4
    ijk(1, 3) = 5
    ijk(3, 1) = 5
    ijk(2, 3) = 6
    ijk(3, 2) = 6
!
! -- On convertit le tenseur
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    ftos(ijk(i, j), ijk(k, l)) = ftot(i, j, k, l)
                end do
            end do
        end do
    end do
!
! -- On rajoute les racines de 2
!
    do j = 4, 6
        do i = 1, 6
            ftos(i, j) = ftos(i, j)*rac2
        end do
    end do
!
    do j = 1, 6
        do i = 4, 6
            ftos(i, j) = ftos(i, j)*rac2
        end do
    end do
!
end subroutine
