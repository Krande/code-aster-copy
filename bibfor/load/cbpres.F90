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
subroutine cbpres(load, mesh, model, geomDime, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/cafotu.h"
#include "asterfort/capres.h"
!
    character(len=8), intent(in)  :: load, mesh, model
    integer(kind=8), intent(in)           :: geomDime
    character(len=4), intent(in)  :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads PRES_REP / FORCE_TUYAU
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbOccPresRep, nbOccForceTuyau
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
    call getfac(keywordfact, nbOccPresRep)
    if (nbOccPresRep .ne. 0) then
        call capres(load, mesh, model, geomDime, valeType, nbOccPresRep)
        mapAlreadyCreated = ASTER_TRUE
    end if
!
! - FORCE_TUYAU loading
!
    keywordfact = 'FORCE_TUYAU'
    call getfac(keywordfact, nbOccForceTuyau)
    if (nbOccForceTuyau .ne. 0) then
        call cafotu(load, model, mapAlreadyCreated, mesh, geomDime, valeType, nbOccForceTuyau)
    end if
!
end subroutine
