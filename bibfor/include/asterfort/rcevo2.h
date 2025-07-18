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
    subroutine rcevo2(nbinti, kinti, csigm, cinst, csiex,&
                      kemixt, cstex, csmex, lfatig, flexio,&
                      lrocht, cnoc, cresu, cpres, lsymm)
        integer(kind=8) :: nbinti
        character(len=16) :: kinti
        character(len=24) :: csigm
        character(len=24) :: cinst
        character(len=24) :: csiex
        aster_logical :: kemixt
        character(len=24) :: cstex
        character(len=24) :: csmex
        aster_logical :: lfatig
        aster_logical :: flexio
        aster_logical :: lrocht
        character(len=24) :: cnoc
        character(len=24) :: cresu
        character(len=24) :: cpres
        aster_logical :: lsymm
    end subroutine rcevo2
end interface
