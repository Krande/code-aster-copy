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
    subroutine gkmet3(nnoff, chfond, iadrgk, milieu, connex,&
                      iadgks, iadgki, abscur, num, typdis)
        integer(kind=8)           :: nnoff
        character(len=24) :: chfond
        integer(kind=8)           :: iadrgk
        aster_logical     :: milieu
        aster_logical     :: connex
        integer(kind=8)           :: iadgks
        integer(kind=8)           :: iadgki
        character(len=24) :: abscur
        integer(kind=8)           :: num
        character(len=16) :: typdis
    end subroutine gkmet3
end interface
