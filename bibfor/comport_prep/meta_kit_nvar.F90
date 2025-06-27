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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine meta_kit_nvar(metaPhas, metaRela, metaGlob, &
                         nbMetaPhas, nbVariMetaRela, nbVariMetaGlob)
!
    implicit none
!
#include "asterc/lccree.h"
#include "asterc/lcinfo.h"
#include "asterc/lcdiscard.h"
!
    character(len=16), intent(in) :: metaPhas, metaRela, metaGlob
    integer(kind=8), intent(out) :: nbMetaPhas, nbVariMetaRela, nbVariMetaGlob
!
! --------------------------------------------------------------------------------------------------
!
! META_*
!
! Number of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  metaPhas         : name of phase
! In  metaRela         : name of relation for each phase
! In  metaGlob         : name of relation for global behaviour
! Out nbMetaPhas       : number of phases
! Out nbVariMetaRela   : number of internal state variables for each phase
! Out nbVariMetaGlob   : number of internal state variables for global behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: metaPhasPy, metaRelaPy, metaGlobPy
    integer(kind=8) :: idummy, idummy2
!
! --------------------------------------------------------------------------------------------------
!
    nbMetaPhas = 0
    nbVariMetaRela = 0
    nbVariMetaGlob = 0
    call lccree(1, metaPhas, metaPhasPy)
    call lccree(1, metaRela, metaRelaPy)
    call lccree(1, metaGlob, metaGlobPy)
    call lcinfo(metaPhasPy, idummy, nbMetaPhas, idummy2)
    call lcinfo(metaRelaPy, idummy, nbVariMetaRela, idummy2)
    call lcinfo(metaGlobPy, idummy, nbVariMetaGlob, idummy2)
    call lcdiscard(metaPhasPy)
    call lcdiscard(metaRelaPy)
    call lcdiscard(metaGlobPy)
!
end subroutine
