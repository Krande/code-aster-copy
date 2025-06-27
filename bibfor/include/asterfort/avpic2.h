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
    subroutine avpic2(method, nbvec, nbordr, jrtrv, jitrv,&
                      npoin, jvalpo, jvalor, npic, jpic,&
                      jordpi)
        integer(kind=8) :: nbordr
        integer(kind=8) :: nbvec
        character(len=8) :: method
        integer(kind=8) :: jrtrv
        integer(kind=8) :: jitrv
        integer(kind=8) :: npoin(nbvec)
        integer(kind=8) :: jvalpo
        integer(kind=8) :: jvalor
        integer(kind=8) :: npic(nbvec)
        integer(kind=8) :: jpic
        integer(kind=8) :: jordpi
    end subroutine avpic2
end interface
