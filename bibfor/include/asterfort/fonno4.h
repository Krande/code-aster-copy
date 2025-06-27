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
    subroutine fonno4(ndim, macofo, noma, nbmac, tablev,&
                      noe, nbnoff, indic)
        integer(kind=8) :: ndim
        character(len=19) :: macofo
        character(len=8) :: noma
        integer(kind=8) :: nbmac
        integer(kind=8) :: tablev(2)
        integer(kind=8) :: noe(4, 4)
        integer(kind=8) :: nbnoff
        integer(kind=8) :: indic(4)
    end subroutine fonno4
end interface
