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
    subroutine xffcr(nfon, jfono, jbaso, jtailo, jindpt,&
                     typfon, jfon, jnofaf, jbas, jtail)
        integer(kind=8) :: nfon
        integer(kind=8) :: jfono
        integer(kind=8) :: jbaso
        integer(kind=8) :: jtailo
        integer(kind=8) :: jindpt
        character(len=19) :: typfon
        integer(kind=8) :: jfon
        integer(kind=8) :: jnofaf
        integer(kind=8) :: jbas
        integer(kind=8) :: jtail
    end subroutine xffcr
end interface
