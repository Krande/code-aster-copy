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

subroutine pmfmats(nomats)
!
! --------------------------------------------------------------------------------------------------
!
!  Retourne le nom du materiau "section" (MATER_SEC) de l'élement PMF courant
!  Si l'élément n'est pas PMF, retourne : ' '
!
! --------------------------------------------------------------------------------------------------
!   out
!       nomats     : nom du materiau "section"
! --------------------------------------------------------------------------------------------------
!
!
    implicit none

#include "jeveux.h"
#include "MultiFiber_type.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lteatt.h"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=*), intent(out) :: nomats
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: isdcom, nbgfmx
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8), pointer :: cpri(:) => null()
!
! --------------------------------------------------------------------------------------------------
    if (.not. lteatt('TYPMOD2', 'PMF')) then
        nomats = ' '
    else
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jeveuo(compor(MULTCOMP), 'L', isdcom)
        call jeveuo(compor(MULTCOMP) (1:8)//'.CPRI', 'L', vi=cpri)
        ! Nombre de groupe de fibre
        nbgfmx = cpri(MULTI_FIBER_NBGRFIBR)
        ! Le matériau de torsion c'est le dernier
        nomats = zk24(isdcom-1+nbgfmx*MULTI_FIBER_SIZEK+1) (1:8)
    end if

end subroutine
