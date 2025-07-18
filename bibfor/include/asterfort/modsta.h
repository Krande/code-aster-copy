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
    subroutine modsta(motcle, matfac, matpre, solveu, lmatm,&
                      nume, iddl, coef, neq, nbmode,&
                      zrmod)
        integer(kind=8) :: neq
        character(len=*) :: motcle
        character(len=*) :: matfac
        character(len=*) :: matpre
        character(len=*) :: solveu
        integer(kind=8) :: lmatm
        character(len=*) :: nume
        integer(kind=8) :: iddl(*)
        real(kind=8) :: coef(*)
        integer(kind=8) :: nbmode
        real(kind=8) :: zrmod(neq, *)
    end subroutine modsta
end interface
