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
    subroutine preml0(n1, n2, diag, col, delg,&
                      prno, deeq, nec, p, q,&
                      lbd1, lbd2, rl, rl1, rl2,&
                      nrl, lt, lmat)
        integer(kind=8) :: n1
        integer(kind=8) :: n2
        integer(kind=8) :: diag(0:*)
        integer(kind=8) :: col(*)
        integer(kind=8) :: delg(*)
        integer(kind=8) :: prno(*)
        integer(kind=8) :: deeq(*)
        integer(kind=8) :: nec
        integer(kind=8) :: p(*)
        integer(kind=8) :: q(*)
        integer(kind=8) :: lbd1(n1)
        integer(kind=8) :: lbd2(n1)
        integer(kind=8) :: rl(4, *)
        integer(kind=8) :: rl1(*)
        integer(kind=8) :: rl2(*)
        integer(kind=8) :: nrl
        integer(kind=8) :: lt
        integer(kind=8) :: lmat
    end subroutine preml0
end interface
