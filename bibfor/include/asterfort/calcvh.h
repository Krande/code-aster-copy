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
    subroutine calcvh(nbmat, materf, eta, vp, sigeqe,&
                      vh, vg)
        integer(kind=8) :: nbmat
        real(kind=8) :: materf(nbmat, 2)
        real(kind=8) :: eta
        real(kind=8) :: vp(3)
        real(kind=8) :: sigeqe
        real(kind=8) :: vh
        real(kind=8) :: vg
    end subroutine calcvh
end interface
