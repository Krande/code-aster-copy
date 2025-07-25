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
    subroutine mltclj(nb, n, ll, m, it,&
                      p, front, frn, adper, trav,&
                      c)
        integer(kind=8) :: p
        integer(kind=8) :: nb
        integer(kind=8) :: n
        integer(kind=8) :: ll
        integer(kind=8) :: m
        integer(kind=8) :: it
        complex(kind=8) :: front(*)
        complex(kind=8) :: frn(*)
        integer(kind=8) :: adper(*)
        complex(kind=8) :: trav(p, nb, *)
        complex(kind=8) :: c(nb, nb, *)
    end subroutine mltclj
end interface
