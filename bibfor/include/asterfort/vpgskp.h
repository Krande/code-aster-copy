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
    subroutine vpgskp(nbeq, nconv, vect, alpha, lmatb,&
                      typeps, vaux, ddlexc, delta)
        integer(kind=8) :: nconv
        integer(kind=8) :: nbeq
        real(kind=8) :: vect(nbeq, nconv)
        real(kind=8) :: alpha
        integer(kind=8) :: lmatb
        integer(kind=8) :: typeps
        real(kind=8) :: vaux(nbeq)
        integer(kind=8) :: ddlexc(nbeq)
        real(kind=8) :: delta(nconv)
    end subroutine vpgskp
end interface
