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
    subroutine inclis(nomres, ssta, sstb, intfa, intfb,&
                      fmlia, fplian, fplibn, fpliao, fplibo,&
                      iada, iadb, numlis, matprj)
        character(len=8) :: nomres
        character(len=8) :: ssta
        character(len=8) :: sstb
        character(len=8) :: intfa
        character(len=8) :: intfb
        character(len=24) :: fmlia
        character(len=24) :: fplian
        character(len=24) :: fplibn
        character(len=24) :: fpliao
        character(len=24) :: fplibo
        integer(kind=8) :: iada(3)
        integer(kind=8) :: iadb(3)
        integer(kind=8) :: numlis
        character(len=8) :: matprj
    end subroutine inclis
end interface
