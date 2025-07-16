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
subroutine cbpesa(load, mesh, valeType)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/capesa.h"
!
    character(len=8), intent(in) :: load, mesh
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load PESANTEUR
!
! --------------------------------------------------------------------------------------------------
!
! In  load      : load
! In  mesh      : mesh
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
    call getfac('PESANTEUR', nbOcc)
    if (nbOcc .ne. 0) then
        call capesa(load, mesh, valeType, nbOcc)
    end if
end subroutine
