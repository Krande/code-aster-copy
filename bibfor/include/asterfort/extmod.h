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
interface
    subroutine extmod(basemo, numddl, nume, nbnumo, dmode,&
                      nbeq, nbnoe, iddl, nbddl)
        integer(kind=8) :: nbddl
        integer(kind=8) :: nbnoe
        integer(kind=8) :: nbnumo
        character(len=8) :: basemo
        character(len=14) :: numddl
        integer(kind=8) :: nume(nbnumo)
        real(kind=8) :: dmode(nbddl*nbnoe*nbnumo)
        integer(kind=8) :: nbeq
        integer(kind=8) :: iddl(nbddl)
    end subroutine extmod
end interface
