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
    subroutine xpoajn(maxfem, ino, lsn, jdirno, prefno,&
                      nfiss, he, nnn, inn, inntot,&
                      nbnoc, nbnofi, inofi, co, iacoo2)
        integer(kind=8) :: nfiss
        character(len=8) :: maxfem
        integer(kind=8) :: ino
        real(kind=8) :: lsn(nfiss)
        integer(kind=8) :: jdirno
        character(len=2) :: prefno(4)
        integer(kind=8) :: he(nfiss)
        integer(kind=8) :: nnn
        integer(kind=8) :: inn
        integer(kind=8) :: inntot
        integer(kind=8) :: nbnoc
        integer(kind=8) :: nbnofi
        integer(kind=8) :: inofi
        real(kind=8) :: co(3)
        integer(kind=8) :: iacoo2
    end subroutine xpoajn
end interface
