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
    subroutine mdexch(nofimd, idfimd, nochmd, numpt, numord,&
                      nbcmpc, nomcmc, nbvato, typent, typgeo,&
                      existc, nbcmfi, nmcmfi, nbval, nbprof,&
                      codret)
        character(len=*) :: nofimd
        med_idt :: idfimd
        character(len=*) :: nochmd
        integer(kind=8) :: numpt
        integer(kind=8) :: numord
        integer(kind=8) :: nbcmpc
        character(len=*) :: nomcmc
        integer(kind=8) :: nbvato
        integer(kind=8) :: typent
        integer(kind=8) :: typgeo
        integer(kind=8) :: existc
        integer(kind=8) :: nbcmfi
        character(len=*) :: nmcmfi
        integer(kind=8) :: nbval
        integer(kind=8) :: nbprof
        integer(kind=8) :: codret
    end subroutine mdexch
end interface
