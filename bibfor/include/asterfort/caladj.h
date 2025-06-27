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
    subroutine caladj(col, diag, xadj, adjncy, n,&
                      nnz, deb, tab, suiv, lmat,&
                      ladjn, nrl)
        integer(kind=8) :: lmat
        integer(kind=8) :: n
        integer(kind=8) :: col(lmat)
        integer(kind=8) :: diag(0:n)
        integer(kind=8) :: xadj(n+1)
        integer(kind=8) :: adjncy(*)
        integer(kind=8) :: nnz(n)
        integer(kind=8) :: deb(n)
        integer(kind=8) :: tab(*)
        integer(kind=8) :: suiv(*)
        integer(kind=8) :: ladjn
        integer(kind=8) :: nrl
    end subroutine caladj
end interface
