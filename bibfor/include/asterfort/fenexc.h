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
    subroutine fenexc(noma, nomnoa, long, nbn, nuno,&
                      diax, nbnfen, noefen, disfen)
        character(len=8) :: noma
        character(len=8) :: nomnoa
        real(kind=8) :: long
        integer(kind=8) :: nbn
        integer(kind=8) :: nuno(*)
        real(kind=8) :: diax(*)
        integer(kind=8) :: nbnfen
        integer(kind=8) :: noefen(*)
        real(kind=8) :: disfen(*)
    end subroutine fenexc
end interface
