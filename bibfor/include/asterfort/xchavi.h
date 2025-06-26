! --------------------------------------------------------------------
! Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org
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
    subroutine xchavi(actpoi, jbasc, jffis, jfon, jvit,&
                      jbeta, ndim, nfonn, sifval)
        integer(kind=8) :: actpoi
        integer(kind=8) :: jbasc
        integer(kind=8) :: jffis
        integer(kind=8) :: jfon
        integer(kind=8) :: jvit
        integer(kind=8) :: jbeta
        integer(kind=8) :: ndim
        integer(kind=8) :: nfonn
        integer(kind=8) :: sifval
    end subroutine xchavi
end interface
