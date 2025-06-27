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
    subroutine aceat2(nbtuy, eltuy, notuy, nbpart, noex1,&
                      noex2, nbmap, elpar, nopar, nno)
        integer(kind=8) :: nno
        integer(kind=8) :: nbpart
        integer(kind=8) :: nbtuy
        integer(kind=8) :: eltuy(nbtuy)
        integer(kind=8) :: notuy(nno*nbtuy)
        integer(kind=8) :: noex1(nbpart)
        integer(kind=8) :: noex2(nbpart)
        integer(kind=8) :: nbmap(nbpart)
        integer(kind=8) :: elpar(nbpart, nbtuy)
        integer(kind=8) :: nopar(nbpart, nno, nbtuy)
    end subroutine aceat2
end interface
