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
    subroutine calcfe(nr, ndt, nvi, vind, df,&
                      gamsns, fe, fp, iret)
        integer(kind=8) :: nr
        integer(kind=8) :: ndt
        integer(kind=8) :: nvi
        real(kind=8) :: vind(*)
        real(kind=8) :: df(3, 3)
        real(kind=8) :: gamsns(3, 3)
        real(kind=8) :: fe(3, 3)
        real(kind=8) :: fp(3, 3)
        integer(kind=8) :: iret
    end subroutine calcfe
end interface
