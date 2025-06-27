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
#include "asterf_types.h"
!
interface 
    subroutine cremnl(reprise, baseno, numrep, nbordr0, nbordr, nbpt, neq,&
                      nbhar, imat, numedd, parcho, nbchoc, vk8, modrep)
        integer(kind=8) :: nbhar
        integer(kind=8) :: numrep
        integer(kind=8) :: neq
        aster_logical :: reprise
        character(len=8) :: baseno
        integer(kind=8) :: nbordr0
        integer(kind=8) :: nbordr
        integer(kind=8) :: nbpt
        integer(kind=8) :: imat(2)
        character(len=24) :: numedd
        character(len=14) :: parcho
        integer(kind=8) :: nbchoc
        character(len=8) :: vk8
        character(len=8) :: modrep
    end subroutine cremnl
end interface 
