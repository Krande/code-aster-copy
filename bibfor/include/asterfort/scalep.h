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
    subroutine scalep(spectr, noma, base, nuor, nbm,&
                      imodi, nbmr, nbexcp, ltable, iaxe,&
                      scal)
        integer(kind=8) :: nbexcp
        integer(kind=8) :: nbmr
        integer(kind=8) :: nbm
        character(len=19) :: spectr
        character(len=8) :: noma
        character(len=19) :: base
        integer(kind=8) :: nuor(nbm)
        integer(kind=8) :: imodi
        aster_logical :: ltable
        integer(kind=8) :: iaxe
        real(kind=8) :: scal(nbexcp, nbmr)
    end subroutine scalep
end interface
