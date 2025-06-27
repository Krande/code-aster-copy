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
!
interface
    subroutine finlf2(ndim, delta, alpha, kn, kt,&
                      mu, Bn, Bt, m1, m2, cbar, d1, res)
        integer(kind=8) :: ndim
        real(kind=8) :: delta(ndim)
        real(kind=8) :: alpha
        real(kind=8) :: kn
        real(kind=8) :: kt
        real(kind=8) :: mu
        real(kind=8) :: Bn
        real(kind=8) :: Bt
        real(kind=8) :: m1
        real(kind=8) :: m2
        real(kind=8) :: cbar
        real(kind=8) :: d1
        real(kind=8) :: res
    end subroutine finlf2
end interface
