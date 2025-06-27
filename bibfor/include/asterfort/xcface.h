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
    subroutine xcface(lsn, lst, jgrlsn, igeom,&
                      enr, nfiss, ifiss, fisco, nfisc,&
                      noma, nmaabs, typdis, pinter, ninter, ainter,&
                      nface, nptf, cface, minlst)
        integer(kind=8) :: nfisc
        character(len=8) :: elref
        real(kind=8) :: lsn(*)
        real(kind=8) :: lst(*)
        integer(kind=8) :: jgrlsn
        integer(kind=8) :: igeom
        character(len=16) :: enr
        integer(kind=8) :: nfiss
        integer(kind=8) :: ifiss
        integer(kind=8) :: fisco(*)
        character(len=8) :: noma
        integer(kind=8) :: nmaabs
        character(len=16) :: typdis
        real(kind=8) :: pinter(*)
        integer(kind=8) :: ninter
        real(kind=8) :: ainter(*)
        integer(kind=8) :: nface
        integer(kind=8) :: nptf
        integer(kind=8) :: cface(30, 6)
        real(kind=8) :: minlst
    end subroutine xcface
end interface
