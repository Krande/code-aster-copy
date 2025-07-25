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
    subroutine fonbas2(noma, basnof, typm, fonoeu, coorfond, nbnoff, absfon,&
                      basloc, abscur, lnno, ltno)
        character(len=8)  :: noma
        character(len=19) :: basnof
        character(len=8)  :: typm
        character(len=24) :: fonoeu
        character(len=24) :: coorfond
        integer(kind=8)           :: nbnoff
        character(len=24) :: absfon
        character(len=19) :: basloc
        character(len=24) :: abscur
        character(len=19) :: lnno
        character(len=19) :: ltno
    end subroutine fonbas2
end interface
