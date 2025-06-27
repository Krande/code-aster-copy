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
    subroutine mnlqd2(ind, imat, neq, ninc, nd,&
                      nchoc, h, hf, parcho, xcdl,&
                      adime, xvect, xtemp)
        integer(kind=8) :: hf
        integer(kind=8) :: ninc
        integer(kind=8) :: neq
        integer(kind=8) :: ind
        integer(kind=8) :: imat(2)
        integer(kind=8) :: nd
        integer(kind=8) :: nchoc
        integer(kind=8) :: h
        character(len=14) :: parcho
        character(len=14) :: xcdl
        character(len=14) :: adime
        character(len=14) :: xvect
        character(len=14) :: xtemp
    end subroutine mnlqd2
end interface 
