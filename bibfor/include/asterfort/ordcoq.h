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
    subroutine ordcoq(imod, nbm, icoq, nbno, numno,&
                      inomax, nbnoto, coordo, iaxe, defm,&
                      nunoe0, drmax, torco)
        integer(kind=8) :: nbnoto
        integer(kind=8) :: nbno
        integer(kind=8) :: nbm
        integer(kind=8) :: imod
        integer(kind=8) :: icoq
        integer(kind=8) :: numno(nbno)
        integer(kind=8) :: inomax
        real(kind=8) :: coordo(3, nbnoto)
        integer(kind=8) :: iaxe
        real(kind=8) :: defm(2, nbnoto, nbm)
        integer(kind=8) :: nunoe0
        real(kind=8) :: drmax
        real(kind=8) :: torco(4, nbm)
    end subroutine ordcoq
end interface
