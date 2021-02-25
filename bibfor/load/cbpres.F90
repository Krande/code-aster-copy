! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine cbpres(load, mesh, ligrmo, ndim, valeType)
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/cafotu.h"
#include "asterfort/capres.h"
!
character(len=8), intent(in) :: load, mesh
character(len=19), intent(in) :: ligrmo
integer, intent(in) :: ndim
character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads PRES_REP / FORCE_TUYAU
!
! --------------------------------------------------------------------------------------------------
!
!
! In  load      : load
! In  mesh      : mesh
! In  ligrmo    : model <LIGREL>
! In  ndim      : dimension of space
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbocc
    aster_logical :: mapAlreadyCreated
    character(len=16) :: keywordfact
!
! --------------------------------------------------------------------------------------------------
!
    mapAlreadyCreated = ASTER_FALSE
!
! - PRES_REP loading
!
    keywordfact = 'PRES_REP'
    call getfac(keywordfact, nbocc)
    if (nbocc .ne. 0) then
        call capres(load, ligrmo, mesh, ndim, valeType, nbocc)
        mapAlreadyCreated = ASTER_TRUE
    endif
!
! - FORCE_TUYAU loading
!
    keywordfact = 'FORCE_TUYAU'
    call getfac(keywordfact, nbocc)
    if (nbocc .ne. 0) then
        call cafotu(load, ligrmo, mapAlreadyCreated, mesh, ndim, valeType, nbocc)
    endif
!
end subroutine
