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
    subroutine lrcnme(chanom, nochmd, nomamd, nomaas, nomgd,&
                      typent, nbcmpv, ncmpva, ncmpvm, iinst,&
                      numpt, numord, inst, crit, prec,&
                      nrofic, codret, base)
        character(len=*) :: chanom
        character(len=*) :: nochmd
        character(len=*) :: nomamd
        character(len=8) :: nomaas
        character(len=8) :: nomgd
        integer(kind=8) :: typent
        integer(kind=8) :: nbcmpv
        character(len=*) :: ncmpva
        character(len=*) :: ncmpvm
        integer(kind=8) :: iinst
        integer(kind=8) :: numpt
        integer(kind=8) :: numord
        real(kind=8) :: inst
        character(len=8) :: crit
        real(kind=8) :: prec
        integer(kind=8) :: nrofic
        integer(kind=8) :: codret
        character(len=1), optional, intent(in) :: base
    end subroutine lrcnme
end interface
