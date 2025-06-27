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
    subroutine cordd2(jprn1, jprn2, ili, ecodl, nec,&
                      ncmp, n, nddloc, pos)
        integer(kind=8) :: jprn1
        integer(kind=8) :: jprn2
        integer(kind=8) :: ili
        integer(kind=8) :: ecodl(*)
        integer(kind=8) :: nec
        integer(kind=8) :: ncmp
        integer(kind=8) :: n
        integer(kind=8) :: nddloc
        integer(kind=8) :: pos(1)
    end subroutine cordd2
end interface
