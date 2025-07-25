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
    subroutine pacoa2(lisi1z, lisi2z, lonli1, lonli2, noma1z,&
                      noma2z, liso1z, liso2z, lonlis)
        character(len=*) :: lisi1z
        character(len=*) :: lisi2z
        integer(kind=8) :: lonli1
        integer(kind=8) :: lonli2
        character(len=*) :: noma1z
        character(len=*) :: noma2z
        character(len=*) :: liso1z
        character(len=*) :: liso2z
        integer(kind=8) :: lonlis
    end subroutine pacoa2
end interface
