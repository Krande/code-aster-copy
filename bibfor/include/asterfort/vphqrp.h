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
    subroutine vphqrp(mat, neq, mxeq, icode, w,&
                      z, iz, wk, mxiter, ier,&
                      nitqr)
        integer(kind=8) :: mxeq
        integer(kind=8) :: neq
        real(kind=8) :: mat(mxeq, 1)
        integer(kind=8) :: icode
        real(kind=8) :: w(1)
        real(kind=8) :: z(1)
        integer(kind=8) :: iz
        real(kind=8) :: wk(neq, 1)
        integer(kind=8) :: mxiter
        integer(kind=8) :: ier
        integer(kind=8) :: nitqr
    end subroutine vphqrp
end interface
