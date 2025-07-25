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
    subroutine ccl21j(fronti, frontj, frn, j, l,&
                      n, n1, t1, t2)
        integer(kind=8) :: n
        complex(kind=8) :: fronti(*)
        complex(kind=8) :: frontj(*)
        complex(kind=8) :: frn(*)
        integer(kind=8) :: j
        integer(kind=8) :: l
        integer(kind=8) :: n1
        complex(kind=8) :: t1(n)
        complex(kind=8) :: t2(n)
    end subroutine ccl21j
end interface
