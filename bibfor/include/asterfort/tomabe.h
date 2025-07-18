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
    subroutine tomabe(chmat, nmabet, nbmabe, mailla, nbnoma,&
                      mail2d, nbnobe, nunobe, xflu, xret,&
                      regl)
        character(len=8) :: chmat
        character(len=24) :: nmabet
        integer(kind=8) :: nbmabe
        character(len=8) :: mailla
        integer(kind=8) :: nbnoma
        aster_logical :: mail2d
        integer(kind=8) :: nbnobe
        character(len=19) :: nunobe
        real(kind=8) :: xflu
        real(kind=8) :: xret
        character(len=4) :: regl
    end subroutine tomabe
end interface
