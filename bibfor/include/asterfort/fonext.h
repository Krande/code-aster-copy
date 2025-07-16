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
    subroutine fonext(noma, cnxinv, jbasno, inoext, inoseg,&
                      nbnoff, jborl, jdirol, jnvdir, iseg)
        character(len=8) :: noma
        character(len=19) :: cnxinv
        integer(kind=8) :: jbasno
        integer(kind=8) :: inoext
        integer(kind=8) :: inoseg
        integer(kind=8) :: nbnoff
        integer(kind=8) :: jborl
        integer(kind=8) :: jdirol
        integer(kind=8) :: jnvdir
        integer(kind=8) :: iseg
    end subroutine fonext
end interface
