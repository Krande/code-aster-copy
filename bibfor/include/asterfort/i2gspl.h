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
    subroutine i2gspl(debspl, tvois1, tvois2, tplace, schm,&
                      achm, pts, pta)
        integer(kind=8) :: debspl
        integer(kind=8) :: tvois1(*)
        integer(kind=8) :: tvois2(*)
        aster_logical :: tplace(*)
        integer(kind=8) :: schm(*)
        integer(kind=8) :: achm(*)
        integer(kind=8) :: pts
        integer(kind=8) :: pta
    end subroutine i2gspl
end interface
