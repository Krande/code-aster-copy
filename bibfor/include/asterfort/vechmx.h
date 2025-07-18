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
    subroutine vechmx(nomo, lischa, ichar, nbch, nomlis,&
                      nbin_maxi, lpain, lchin, lastin, vecele)
        integer(kind=8) :: nbin_maxi
        character(len=8) :: nomo
        character(len=19) :: lischa
        integer(kind=8) :: ichar
        integer(kind=8) :: nbch
        character(len=24) :: nomlis
        character(len=8) :: lpain(nbin_maxi)
        character(len=19) :: lchin(nbin_maxi)
        integer(kind=8) :: lastin
        character(len=19) :: vecele
    end subroutine vechmx
end interface
