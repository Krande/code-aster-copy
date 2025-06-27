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
    subroutine aceat3(noma, nomu, nbtuy, nbpart, nbmap,&
                      elpar, nopar, ivr, nbzk,&
                      nozk, cozk, isens, coor, epsi,&
                      crit, nno, nmmt)
        integer(kind=8) :: nno
        integer(kind=8) :: nbzk
        integer(kind=8) :: nbpart
        integer(kind=8) :: nbtuy
        character(len=8) :: noma
        character(len=8) :: nomu
        integer(kind=8) :: nbmap(nbpart)
        integer(kind=8) :: elpar(nbpart, nbtuy)
        integer(kind=8) :: nopar(nbpart, nno, nbtuy)
        integer(kind=8) :: ivr(*)
        integer(kind=8) :: ifm
        integer(kind=8) :: nozk(nbzk)
        real(kind=8) :: cozk(3*nbzk)
        integer(kind=8) :: isens(nbpart)
        real(kind=8) :: coor(*)
        real(kind=8) :: epsi
        character(len=8) :: crit
        integer(kind=8) :: nmmt(*)
    end subroutine aceat3
end interface
