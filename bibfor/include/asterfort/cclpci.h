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
    subroutine cclpci(option, modele, resuin, resuou, mater , mateco, &
                      carael, ligrel, numord, nbpain, lipain,&
                      lichin, codret)
        character(len=16) :: option
        character(len=8) :: modele
        character(len=8) :: resuin
        character(len=8) :: resuou
        character(len=8) :: mater, mateco
        character(len=8) :: carael
        character(len=24) :: ligrel
        integer(kind=8) :: numord
        integer(kind=8) :: nbpain
        character(len=8) :: lipain(*)
        character(len=24) :: lichin(*)
        integer(kind=8) :: codret
    end subroutine cclpci
end interface
