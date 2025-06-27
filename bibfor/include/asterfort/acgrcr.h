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
    subroutine acgrcr(nbvec, jvectn, jvectu, jvectv, nbordr, &
                      kwork, sompgw, jrwork, tspaq, ipg, &
                      nommet,jvecno, jnorma, forcri,nompar,&
                      vanocr, respc,vnmax)
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
        character(len=16) :: nommet
        integer(kind=8) :: jvecno
        integer(kind=8) :: jnorma
        character(len=16) :: forcri
        character(len=8) ::nompar(35)
        real(kind=8) ::vanocr(23)
        real(kind=8) :: respc(24)
        real(kind=8) :: vnmax(6)
    end subroutine acgrcr
end interface
