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
    subroutine vdrep2(alpha, beta, zilzi, zrlzr, matevn,&
                      matevg)
        real(kind=8) :: alpha
        real(kind=8) :: beta
        integer(kind=8) :: zilzi(*)
        real(kind=8) :: zrlzr(*)
        real(kind=8) :: matevn(2, 2, 1)
        real(kind=8) :: matevg(2, 2, 1)
    end subroutine vdrep2
end interface
