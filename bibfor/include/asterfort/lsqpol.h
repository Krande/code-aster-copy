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
    subroutine lsqpol(ordre, e1, npt, xx, yy,&
                      ordok, poly, sigma)
        integer(kind=8) :: npt
        integer(kind=8) :: ordre
        real(kind=8) :: e1
        real(kind=8) :: xx(npt)
        real(kind=8) :: yy(npt)
        integer(kind=8) :: ordok
        real(kind=8) :: poly(ordre+1)
        real(kind=8) :: sigma
    end subroutine lsqpol
end interface
