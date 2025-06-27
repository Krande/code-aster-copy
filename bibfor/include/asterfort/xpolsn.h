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
    subroutine xpolsn(elrefp, ino, n, jlsn, jlst,&
                      ima, iad, igeom, nfiss, ndime,&
                      ndim, jconx1, jconx2, fisco, co,&
                      lsn, lst)
        integer(kind=8) :: nfiss
        integer(kind=8) :: n
        character(len=8) :: elrefp
        integer(kind=8) :: ino
        integer(kind=8) :: jlsn
        integer(kind=8) :: jlst
        integer(kind=8) :: ima
        integer(kind=8) :: iad
        integer(kind=8) :: igeom
        integer(kind=8) :: ndime
        integer(kind=8) :: ndim
        integer(kind=8) :: jconx1
        integer(kind=8) :: jconx2
        integer(kind=8) :: fisco(*)
        real(kind=8) :: co(3)
        real(kind=8) :: lsn(nfiss)
        real(kind=8) :: lst(nfiss)
    end subroutine xpolsn
end interface
