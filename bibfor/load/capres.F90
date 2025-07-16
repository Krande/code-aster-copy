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
subroutine capres(load, mesh, model, geomDime, valeType, nbOccPresRep)
!
    implicit none
!
#include "asterfort/capres_skin.h"
#include "asterfort/capres_volu.h"
#include "asterfort/dismoi.h"
!
    character(len=8), intent(in)  :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in)  :: valeType
    integer(kind=8), intent(in) :: nbOccPresRep
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads PRES_REP
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
! In  nbOccPresRep     : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: answer
!
! --------------------------------------------------------------------------------------------------
!

! - For skin elements
    call capres_skin(load, mesh, model, geomDime, valeType, nbOccPresRep)

! - For volumic elements
    call dismoi('EXI_COQSOL', model, 'MODELE', repk=answer)
    if (answer .eq. 'OUI') then
        call capres_volu(load, mesh, valeType, nbOccPresRep)
    end if
!
end subroutine
