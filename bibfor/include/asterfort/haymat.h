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
    subroutine haymat(fami, kpg, ksp, mod, imat,&
                      nmat, poum, coefel, coefpl, nvi,&
                      nr)
        integer(kind=8) :: nmat
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=8) :: mod
        integer(kind=8) :: imat
        character(len=*) :: poum
        real(kind=8) :: coefel(nmat)
        real(kind=8) :: coefpl(nmat)
        integer(kind=8) :: nvi
        integer(kind=8) :: nr
    end subroutine haymat
end interface
