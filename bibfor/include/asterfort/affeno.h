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
    subroutine affeno(ioc, ino, nocmp, nbcmp, ncmpgd,&
                      ncmpmx, val, kval, desc, valglo,&
                      kvalgl, type, nec)
        integer(kind=8) :: ioc
        integer(kind=8) :: ino
        character(len=8) :: nocmp(*)
        integer(kind=8) :: nbcmp
        character(len=8) :: ncmpgd(*)
        integer(kind=8) :: ncmpmx
        real(kind=8) :: val(*)
        character(len=8) :: kval(*)
        integer(kind=8) :: desc(*)
        real(kind=8) :: valglo(*)
        character(len=8) :: kvalgl(*)
        character(len=*) :: type
        integer(kind=8) :: nec
    end subroutine affeno
end interface
