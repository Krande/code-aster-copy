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
    subroutine cakg2d(optioz, result, modele, depla, theta,&
                      mate, mateco, lischa, symech, fondf, noeud, &
                      time, iord, nbprup, noprup, &
                      lmoda, puls, compor)
        character(len=16) :: optioz
        character(len=8) :: result
        character(len=8) :: modele
        character(len=24) :: depla
        character(len=24) :: theta
        character(len=24) :: mate, mateco
        character(len=19) :: lischa
        character(len=8) :: symech
        character(len=8) :: fondf
        character(len=8) :: noeud
        real(kind=8) :: time
        integer(kind=8) :: iord
        integer(kind=8) :: nbprup
        character(len=16) :: noprup(*)
        aster_logical :: lmoda
        real(kind=8) :: puls
        character(len=24) :: compor
    end subroutine cakg2d
end interface
