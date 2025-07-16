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
! Person in charge: mickael.abbas at edf.fr
!
subroutine calirc(phenomZ, load, model)
!
    use LoadKinematic_module
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
!
    character(len=*), intent(in) :: phenomZ
    character(len=8), intent(in) :: load, model
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'LIAISON_MAIL'
!
! --------------------------------------------------------------------------------------------------
!
! In  phenom           : phenomenon (MECANIQUE/THERMIQUE/ACOUSTIQUE)
! In  load             : load
! In  model            : model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: valeType = 'REEL'
    character(len=16), parameter :: factorKeyword = 'LIAISON_MAIL'
    character(len=19) :: listLineRela
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbOcc
    aster_logical :: lVerbose
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    lVerbose = (niv .ge. 2)

    call getfac(factorKeyword, nbOcc)
    if (nbOcc .ne. 0) then
        if (phenomZ .eq. 'MECANIQUE') then
            call kineLoadGlueMeshMeca(load, model, valeType, lVerbose, listLineRela)
        elseif (phenomZ .eq. 'THERMIQUE') then
            call kineLoadGlueMeshTher(model, valeType, listLineRela)
        else
            ASSERT(ASTER_FALSE)
        end if
        call aflrch(listLineRela, load, 'LIN', detr_lisrez=ASTER_TRUE)
    end if
!
end subroutine
