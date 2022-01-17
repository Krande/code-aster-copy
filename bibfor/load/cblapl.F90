! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine cblapl(load, mesh, model)
!
implicit none
!
#include "asterc/getfac.h"
#include "asterfort/calapl.h"
!
character(len=8), intent(in) :: load, mesh, model
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load INTE_ELEC
!
! --------------------------------------------------------------------------------------------------
!
! In  load      : load
! In  mesh      : mesh
! In  model     : model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'INTE_ELEC'
    integer :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
    call getfac(keywordfact, nbOcc)
    if (nbOcc .ne. 0) then
        call calapl(load, mesh, model, nbOcc)
    endif
!
end subroutine
