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
    subroutine vpini1(eigsol, modes, solveu, typcon, vecblo, veclag, vecrig,&
                      matpsc, matopa, iretr, nblagr, neqact, npivot, nstoc, omemax, omemin, omeshi,&
                      sigma, mod45, nfreq_calibr_)
        character(len=19) , intent(in)  :: eigsol
        character(len=8)  , intent(in)  :: modes
        character(len=19) , intent(in)  :: solveu
        character(len=16) , intent(in)  :: typcon
        character(len=24) , intent(in)  :: vecblo
        character(len=24) , intent(in)  :: veclag
        character(len=24) , intent(in)  :: vecrig
!!
        character(len=19) , intent(inout)  :: matpsc
        character(len=19) , intent(inout)  :: matopa
        integer(kind=8)           , intent(out) :: iretr
        integer(kind=8), optional , intent(out) :: nfreq_calibr_
        integer(kind=8)           , intent(out) :: nblagr
        integer(kind=8)           , intent(out) :: neqact
        integer(kind=8)           , intent(out) :: npivot
        integer(kind=8)           , intent(out) :: nstoc
        real(kind=8)      , intent(out) :: omemax
        real(kind=8)      , intent(out) :: omemin
        real(kind=8)      , intent(out) :: omeshi
        complex(kind=8)   , intent(out) :: sigma
        character(len=4)  , intent(in)  :: mod45
    end subroutine vpini1
end interface
