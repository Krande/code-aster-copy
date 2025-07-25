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
subroutine kpouli(e, a, nx, l0, l1, &
                  l2, norml1, norml2, amat)
! ......................................................................
!    - FONCTION REALISEE:  CALCUL MATRICE DE RIGIDITE MEPOULI
    implicit none
!                          OPTION : 'FULL_MECA        '
!                          OPTION : 'RIGI_MECA_TANG   '
!
!    - ARGUMENTS:
!        DONNEES:
!
! ......................................................................
!
    real(kind=8) :: e, a, nx, l1(3), l2(3), c123(3), c456(3)
    real(kind=8) :: norml1, norml2, l0, coef1, coef2, coef3, coef4, coef5
    real(kind=8) :: amat(*)
    integer(kind=8) :: i, j, imat
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    coef1 = (e*a/l0-nx/norml1)/norml1**2
    coef2 = e*a/(l0*norml1*norml2)
    coef3 = (e*a/l0-nx/norml2)/norml2**2
    coef4 = nx/norml1
    coef5 = nx/norml2
    imat = 0
!
!*** LIGNES 1, 2, 3
!
    do i = 1, 3
        do j = 1, i
            imat = imat+1
            amat(imat) = coef1*l1(i)*l1(j)
        end do
        amat(imat) = amat(imat)+coef4
    end do
!
!*** LIGNES 4, 5, 6
!
    do i = 1, 3
        do j = 1, 3
            imat = imat+1
            amat(imat) = coef2*l2(i)*l1(j)
        end do
        do j = 1, i
            imat = imat+1
            amat(imat) = coef3*l2(i)*l2(j)
        end do
        amat(imat) = amat(imat)+coef5
    end do
!
!*** LIGNES 7, 8, 9
!
    do i = 1, 3
        do j = 1, 3
            imat = imat+1
            amat(imat) = -coef1*l1(i)*l1(j)-coef2*l2(i)*l1(j)
            if (j .eq. i) then
                amat(imat) = amat(imat)-coef4
            end if
            c123(j) = amat(imat)
        end do
        do j = 1, 3
            imat = imat+1
            amat(imat) = -coef2*l1(i)*l2(j)-coef3*l2(i)*l2(j)
            if (j .eq. i) then
                amat(imat) = amat(imat)-coef5
            end if
            c456(j) = amat(imat)
        end do
        do j = 1, i
            imat = imat+1
            amat(imat) = -c123(j)-c456(j)
        end do
    end do
end subroutine
