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
    subroutine concom(macor, nblir, macoc, nblic, nbnoco,&
                      nococ)
        integer(kind=8) :: nblic
        integer(kind=8) :: nblir
        character(len=8) :: macor(nblir+2)
        character(len=8) :: macoc(nblic+2)
        integer(kind=8) :: nbnoco
        integer(kind=8) :: nococ(*)
    end subroutine concom
end interface
