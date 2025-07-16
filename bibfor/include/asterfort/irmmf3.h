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
interface
    subroutine irmmf3(fid, nomamd, typent, nbrent, nbgrou,&
                      nomgen, nbec, nomast, prefix, typgeo,&
                      nomtyp, nmatyp, nufaen, nufacr, nogrfa,&
                      nofaex, tabaux, infmed, ifm, nosdfu)
        integer(kind=8) :: nbgrou
        integer(kind=8) :: nbrent
        med_idt :: fid
        character(len=*) :: nomamd
        integer(kind=8) :: typent
        character(len=24) :: nomgen(*)
        integer(kind=8) :: nbec
        character(len=8) :: nomast
        character(len=6) :: prefix
        integer(kind=8) :: typgeo(*)
        character(len=8) :: nomtyp(*)
        integer(kind=8) :: nmatyp(*)
        integer(kind=8) :: nufaen(nbrent)
        integer(kind=8) :: nufacr(nbrent)
        character(len=80) :: nogrfa(nbgrou)
        character(len=*) :: nofaex(*)
        integer(kind=8) :: tabaux(*)
        integer(kind=8) :: infmed
        integer(kind=8) :: ifm
        character(len=8) :: nosdfu
    end subroutine irmmf3
end interface
