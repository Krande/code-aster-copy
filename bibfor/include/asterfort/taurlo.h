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
    subroutine taurlo(nbvec, jvectn, jvectu, jvectv, nbordr,&
                      kwork, sompgw, jrwork, tspaq, ipg, dectau, &
                      jvecpg, jnorma)
        integer(kind=8) :: nbvec
        integer(kind=8) :: jvectn
        integer(kind=8) :: jvectu
        integer(kind=8) :: jvectv
        integer(kind=8) :: nbordr
        integer(kind=8) :: kwork
        integer(kind=8) :: sompgw
        integer(kind=8) :: jrwork
        integer(kind=8) :: tspaq
        integer(kind=8) :: ipg
        integer(kind=8) ::dectau
        integer(kind=8) :: jvecpg
        integer(kind=8) :: jnorma
    end subroutine taurlo
end interface
