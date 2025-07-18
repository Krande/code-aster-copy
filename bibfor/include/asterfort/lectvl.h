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
    subroutine lectvl(zcmplx, itype, nbabs, inatur, ideas,&
                      nbmesu, labs, amin, apas, lvalc,&
                      lvalr)
        aster_logical :: zcmplx
        integer(kind=8) :: itype
        integer(kind=8) :: nbabs
        integer(kind=8) :: inatur
        integer(kind=8) :: ideas
        integer(kind=8) :: nbmesu
        integer(kind=8) :: labs
        real(kind=8) :: amin
        real(kind=8) :: apas
        integer(kind=8) :: lvalc
        integer(kind=8) :: lvalr
    end subroutine lectvl
end interface
