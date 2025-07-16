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
interface
    subroutine xfacxt(elp, jpint, jmilt, jnit, jcnset, pinter,&
                      ninter, jphe, ndim, ainter,nface,nptf, cface,&
                      igeom, jlsn, jlst, jaint, jgrlsn)
        integer(kind=8) :: ninter
        integer(kind=8) :: nface
        integer(kind=8) :: cface(30,6)
        integer(kind=8) :: jcnset
        integer(kind=8) :: jnit
        integer(kind=8) :: jmilt
        integer(kind=8) :: jpint
        integer(kind=8) :: nptf
        integer(kind=8) :: ndim
        integer(kind=8) :: jphe
        integer(kind=8) :: igeom
        integer(kind=8) :: jlsn
        integer(kind=8) :: jaint
        integer(kind=8) :: jlst
        integer(kind=8) :: jgrlsn
        real(kind=8) :: pinter(*)
        real(kind=8) :: ainter(*)
        character(len=8) :: elp
    end subroutine xfacxt
end interface
