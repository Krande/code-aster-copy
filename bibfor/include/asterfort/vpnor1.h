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
    subroutine vpnor1(norm, neq, nbmode, ddlexc, vecpro,&
                      isign, numddl, coef)
        integer(kind=8) :: neq
        character(len=*) :: norm
        integer(kind=8) :: nbmode
        integer(kind=8) :: ddlexc(*)
        real(kind=8) :: vecpro(neq, *)
        integer(kind=8) :: isign
        integer(kind=8) :: numddl
        real(kind=8) :: coef(*)
    end subroutine vpnor1
end interface
