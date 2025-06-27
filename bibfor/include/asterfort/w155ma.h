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
    subroutine w155ma(numa, nucou, nicou, nangl, nufib,&
                      motfac, jce2d, jce2l, jce2v, jce5d,&
                      jce5l, jce5v, ksp1, ksp2, c1,&
                      c2, iret)
        integer(kind=8) :: numa
        integer(kind=8) :: nucou
        character(len=3) :: nicou
        integer(kind=8) :: nangl
        integer(kind=8) :: nufib
        character(len=16) :: motfac
        integer(kind=8) :: jce2d
        integer(kind=8) :: jce2l
        integer(kind=8) :: jce2v
        integer(kind=8) :: jce5d
        integer(kind=8) :: jce5l
        integer(kind=8) :: jce5v
        integer(kind=8) :: ksp1
        integer(kind=8) :: ksp2
        real(kind=8) :: c1
        real(kind=8) :: c2
        integer(kind=8) :: iret
    end subroutine w155ma
end interface
