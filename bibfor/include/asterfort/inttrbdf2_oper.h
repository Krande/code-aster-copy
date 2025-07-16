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
    subroutine inttrbdf2_oper(nbequ, par, mgen, kgen, cgen, &
                            ktilda, ftild1, ftild2)
        integer(kind=8)     , intent(in)  :: nbequ
        real(kind=8)              :: par(:)
        real(kind=8), pointer  :: mgen(:)
        real(kind=8), pointer  :: kgen(:)
        real(kind=8), pointer  :: cgen(:)
        real(kind=8), pointer :: ktilda(:)
        real(kind=8), pointer :: ftild1(:)
        real(kind=8), pointer :: ftild2(:)
    end subroutine inttrbdf2_oper
end interface
