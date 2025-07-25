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
subroutine dbrMainPodIncr(lReuse, paraPod, baseOut)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dbr_calcpod_q.h"
#include "asterfort/dbr_calcpod_save.h"
#include "asterfort/dbr_calcpod_size.h"
#include "asterfort/dbr_pod_incr.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    aster_logical, intent(in) :: lReuse
    type(ROM_DS_ParaDBR_POD), intent(in) :: paraPod
    type(ROM_DS_Empi), intent(in) :: baseOut
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUIT
!
! Main subroutine to compute base - For incremental POD
!
! --------------------------------------------------------------------------------------------------
!
! In  lReuse          : .true. if reuse
! In  paraPod          : datastructure for parameters (POD)
! In  baseOut          : output base
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8), pointer :: q(:) => null(), v(:) => null(), s(:) => null()
    integer(kind=8) :: nbMode, nbSnap, m, n
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_59')
    end if
!
! - Get size of snapshots matrix
!
    call dbr_calcpod_size(baseOut, paraPod%snap, &
                          m, n)
!
! - Create snapshots matrix Q
!
    call dbr_calcpod_q(paraPod, baseOut, m, n, q)
!
! - Incremental POD method
!
    call dbr_pod_incr(lReuse, baseOut, paraPod, &
                      q, s, v, nbMode, nbSnap)
!
! - Save base
!
    call dbr_calcpod_save(baseOut, nbMode, nbSnap, s, v)
!
! - Clean
!
    AS_DEALLOCATE(vr=q)
    AS_DEALLOCATE(vr=v)
    AS_DEALLOCATE(vr=s)
!
end subroutine
