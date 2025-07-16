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
subroutine rco3d_elem(nomte, ndim, nddl, &
                      typmaco, nnco, typma3d, nn3d)
!
    implicit none
!
#include "asterf_types.h"

    character(len=16), intent(in) :: nomte
    integer(kind=8), intent(out) :: ndim, nddl, nnco, nn3d
    character(len=8), intent(out) :: typmaco, typma3d

    if (nomte(1:7) .eq. 'RACS2T3') then
        ndim = 3
        typmaco = 'SE2'
        nnco = 2
        typma3d = 'TR3'
        nn3d = 3
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS2T6') then
        ndim = 3
        typmaco = 'SE2'
        nnco = 2
        typma3d = 'TR6'
        nn3d = 6
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS2Q4') then
        ndim = 3
        typmaco = 'SE2'
        nnco = 2
        typma3d = 'QU4'
        nn3d = 4
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS2Q8') then
        ndim = 3
        typmaco = 'SE2'
        nnco = 2
        typma3d = 'QU8'
        nn3d = 8
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS3T3') then
        ndim = 3
        typmaco = 'SE3'
        nnco = 3
        typma3d = 'TR3'
        nn3d = 3
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS3T6') then
        ndim = 3
        typmaco = 'SE3'
        nnco = 3
        typma3d = 'TR6'
        nn3d = 6
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS3Q4') then
        ndim = 3
        typmaco = 'SE3'
        nnco = 3
        typma3d = 'QU4'
        nn3d = 4
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

    if (nomte(1:7) .eq. 'RACS3Q8') then
        ndim = 3
        typmaco = 'SE3'
        nnco = 3
        typma3d = 'QU8'
        nn3d = 8
        nddl = nn3d*ndim+nnco*(ndim+3)
    end if

end subroutine
