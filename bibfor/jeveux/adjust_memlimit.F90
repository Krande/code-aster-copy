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

subroutine adjust_memlimit(show)
    use logging_module, only: DEBUG, LOGLEVEL_MEM, is_enabled
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/utgtme.h"
#include "asterfort/utmess.h"
#include "asterfort/utptme.h"
! ----------------------------------------------------------------------
! Adjust the memory limit = memory for jeveux objects + initial vmpeak
! ----------------------------------------------------------------------
    aster_logical, intent(in) :: show

    integer(kind=8), parameter :: n = 5
    integer(kind=8) :: lfic(n), mfic
    common/fenvje/lfic, mfic

    aster_logical :: init, dbg
    integer(kind=8) :: i, iret
    integer(kind=8), parameter :: npar = 6
    real(kind=8), parameter :: jv_syst = 512.0
    real(kind=8) :: rval(npar), added, sizf
    character(len=8) :: k8tab(npar)
! ----------------------------------------------------------------------
    dbg = is_enabled(LOGLEVEL_MEM, DEBUG)
    k8tab(1) = 'LIMIT_JV'
    k8tab(2) = 'MEM_TOTA'
    k8tab(3) = 'MEM_INIT'
    k8tab(4) = 'VMPEAK'
    k8tab(5) = 'VMSIZE'
    k8tab(6) = 'RLQ_MEM'
    call utgtme(npar, k8tab, rval, iret)
    ASSERT(iret .eq. 0)
    init = rval(3) .lt. 1.
    if (dbg) then
        do i = 1, 6
            write (6, *) "DEBUG: ", k8tab(i), rval(i)
        end do
    end if

    if (rval(4) .le. 0 .or. rval(5) .le. 0) then
        call utmess('I', 'JEVEUX1_75')
    end if

    added = rval(4)
    if (added .gt. 0) then
!       minimal memory to read elements catalogs
        if (init) then
            added = added+jv_syst
            call utptme('LIMIT_JV', rval(1)+jv_syst, iret)
            ASSERT(iret .eq. 0)
        end if
        call utptme('MEM_TOTA', rval(1)+added, iret)
        ASSERT(iret .eq. 0)
        call utptme('MEM_INIT', rval(4), iret)
        ASSERT(iret .eq. 0)
        call utptme('RLQ_MEM ', rval(4), iret)
        ASSERT(iret .eq. 0)
    end if

    if (show .or. dbg) then
        call utgtme(npar, k8tab, rval, iret)
        ASSERT(iret .eq. 0)
    end if
    if (dbg) then
        do i = 1, 6
            write (6, *) "DEBUG: ", k8tab(i), rval(i)
        end do
    end if
    if (show) then
        if (rval(4) .gt. 0) then
            call utmess('I', 'SUPERVIS2_22', nr=3, valr=rval)
        else
            call utmess('I', 'SUPERVIS2_29', nr=2, valr=rval)
        end if
        sizf = mfic/(1024*1024.0d0)
        call utmess('I', 'SUPERVIS2_24', sr=sizf)
    end if

end subroutine
