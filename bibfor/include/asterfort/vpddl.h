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
    subroutine vpddl(raide, masse, neq, nblagr, nbcine,&
                     neqact, dlagr, dbloq, ier)
        integer(kind=8) :: neq
        character(len=19) :: raide
        character(len=19) :: masse
        integer(kind=8) :: nblagr
        integer(kind=8) :: nbcine
        integer(kind=8) :: neqact
        integer(kind=8) :: dlagr(neq)
        integer(kind=8) :: dbloq(neq)
        integer(kind=8) :: ier
    end subroutine vpddl
end interface
