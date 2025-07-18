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

function gtstat(istat)
! person_in_charge: mathieu.courtois at edf.fr
!
!
    use parameters_module, only: ST_OK
    implicit none
#include "asterf_debug.h"
#include "asterf_types.h"
    aster_logical :: gtstat
!     ARGUMENT IN
    integer(kind=8) :: istat
!-----------------------------------------------------------------------
!     FONCTION "GeT STATus" : DIT SI LE STATUT ISTAT EST ACTUELLEMENT
!     ACTIVE OU NON.
!
!     LA VALEUR DU STATUT GLOBAL IGLBST EST STOCKE DANS LE COMMON CGLBST
!-----------------------------------------------------------------------
    integer(kind=8) :: iglbst
    common/cglbst/iglbst
!
    if (istat .eq. ST_OK) then
        gtstat = iglbst .eq. ST_OK
    else
        gtstat = iand(istat, iglbst) .eq. istat
    end if
    DEBUG_MPI('get status: in/returned', istat, gtstat)
end function
