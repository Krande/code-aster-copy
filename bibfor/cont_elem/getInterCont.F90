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
subroutine getInterCont(nbPoinInte, poinInteSlav)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/mesh_pairing_type.h"
#include "jeveux.h"
!
    integer(kind=8), intent(out) :: nbPoinInte
    real(kind=8), intent(out) :: poinInteSlav(2, MAX_NB_INTE)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Quadrature
!
! Get intersection points
!
! --------------------------------------------------------------------------------------------------
!
! Out nbPoinInte       : number of intersection points
! Out poinInteSlav     : coordinates of intersection points (in slave parametric space)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ::  jcont, iPoinInte
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PCONFR', 'L', jcont)
    nbPoinInte = int(zr(jcont-1+1))
    ASSERT(nbPoinInte .le. MAX_NB_INTE)
    poinInteSlav = 0.d0
    do iPoinInte = 1, nbPoinInte
        poinInteSlav(1, iPoinInte) = zr(jcont-1+1+iPoinInte)
        poinInteSlav(2, iPoinInte) = zr(jcont-1+9+iPoinInte)
    end do
!
end subroutine
