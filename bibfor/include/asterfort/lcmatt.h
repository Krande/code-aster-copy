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
    subroutine lcmatt(fami, kpg, ksp, mod, imat,&
                      nmat, poum, rela_comp, coefel, coefpl,&
                      typma, ndt, ndi, nr, nvi)
        integer(kind=8) :: nmat
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=8) :: mod
        integer(kind=8) :: imat
        character(len=*) :: poum
        character(len=16) :: rela_comp
        real(kind=8) :: coefel(nmat)
        real(kind=8) :: coefpl(nmat)
        character(len=8) :: typma
        integer(kind=8) :: ndt
        integer(kind=8) :: ndi
        integer(kind=8) :: nr
        integer(kind=8) :: nvi
    end subroutine lcmatt
end interface
