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
    subroutine prmadj(nbnd, neq, n2, adjncy, xadj,&
                      xadjd, liste, q, noeud)
        integer(kind=8) :: n2
        integer(kind=8) :: neq
        integer(kind=8) :: nbnd
        integer(kind=8) :: adjncy(*)
        integer(kind=8) :: xadj(neq+1)
        integer(kind=8) :: xadjd(*)
        integer(kind=8) :: liste(neq)
        integer(kind=8) :: q(n2)
        integer(kind=8) :: noeud(*)
    end subroutine prmadj
end interface
