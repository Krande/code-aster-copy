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
subroutine cbgrai(load, mesh, model, valeType)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/cagrai.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load PRE_GRAD_TEMP
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'PRE_GRAD_TEMP'
    integer(kind=8) :: nbfac
!
! --------------------------------------------------------------------------------------------------
!
    call getfac(keywordFact, nbfac)
!
    if (nbfac .ne. 0) then
        call cagrai(load, mesh, model, valeType)
    end if
!
end subroutine
