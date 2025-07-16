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
    subroutine irmmfa(fid, nomamd, nbnoeu, nbmail, nomast,&
                      nbgrno, nomgno, nbgrma, nomgma, prefix,&
                      typgeo, nomtyp, nmatyp, infmed, nosdfu)
        med_idt :: fid
        character(len=*) :: nomamd
        integer(kind=8) :: nbnoeu
        integer(kind=8) :: nbmail
        character(len=8) :: nomast
        integer(kind=8) :: nbgrno
        character(len=24) :: nomgno(*)
        integer(kind=8) :: nbgrma
        character(len=24) :: nomgma(*)
        character(len=6) :: prefix
        integer(kind=8) :: typgeo(*)
        character(len=8) :: nomtyp(*)
        integer(kind=8) :: nmatyp(*)
        integer(kind=8) :: infmed
        character(len=8) :: nosdfu
    end subroutine irmmfa
end interface
