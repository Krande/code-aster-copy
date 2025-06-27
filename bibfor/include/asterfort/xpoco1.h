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
    subroutine xpoco1(dirma, nbma, dirno, nbno, ma1,&
                      ma2, jnivgr)
        integer(kind=8) :: nbno
        integer(kind=8) :: nbma
        integer(kind=8) :: dirma(nbma)
        integer(kind=8) :: dirno(nbno)
        character(len=8) :: ma1
        character(len=8) :: ma2
        integer(kind=8) :: jnivgr
    end subroutine xpoco1
end interface
