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
    subroutine nmdeca(sddisc, iterat, ievdac, nomlis, instam,&
                      deltat, nbrpas, dtmin, ldcext, durdec,&
                      retdec)
        character(len=19) :: sddisc
        integer(kind=8) :: iterat
        integer(kind=8) :: ievdac
        character(len=24) :: nomlis
        real(kind=8) :: instam
        real(kind=8) :: deltat
        integer(kind=8) :: nbrpas
        real(kind=8) :: dtmin
        aster_logical :: ldcext
        real(kind=8) :: durdec
        integer(kind=8) :: retdec
    end subroutine nmdeca
end interface
