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
    subroutine fondpl(modele, mate, mateco, numedd, neq, chondp,&
                      nchond, vecond, veonde, vaonde, temps,&
                      foonde)
        integer(kind=8) :: nchond
        integer(kind=8) :: neq
        character(len=24) :: modele
        character(len=8) :: mate
        character(len=24) :: mateco
        character(len=24) :: numedd
        character(len=8) :: chondp(nchond)
        character(len=24) :: vecond
        character(len=24) :: veonde
        character(len=24) :: vaonde
        real(kind=8) :: temps
        real(kind=8) :: foonde(neq)
    end subroutine fondpl
end interface
