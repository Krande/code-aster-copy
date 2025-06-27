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
    subroutine vampli(vwork, tdisp, liste, nbt, nbordr,&
                      numini, nbp, tspaq, nomopt, cxsr)
        integer(kind=8) :: nbp
        integer(kind=8) :: tdisp
        real(kind=8) :: vwork(tdisp)
        integer(kind=8) :: liste(nbp)
        integer(kind=8) :: nbt
        integer(kind=8) :: nbordr
        integer(kind=8) :: numini
        integer(kind=8) :: tspaq
        character(len=16) :: nomopt
        character(len=19) :: cxsr
    end subroutine vampli
end interface
