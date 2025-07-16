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

subroutine op0165()
    implicit none
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!     OPERATEUR POST_RCCM
!
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infmaj.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/rccome.h"
#include "asterfort/utmess.h"
#include "asterfort/rcevol.h"
#include "asterfort/rc3600.h"
#include "asterfort/rc3600_momeq.h"
#include "asterfort/rc3200.h"
#include "asterfort/titre.h"
!     ------------------------------------------------------------------
!
    real(kind=8) :: symax
    character(len=16) :: typtab, typmec, kopt(4)
    character(len=24) :: option
    integer(kind=8) :: n1, nbopt, icodre
    character(len=8) :: nommat, symm
    aster_logical :: lsymm
!
! DEB ------------------------------------------------------------------
!
    call infmaj()
!
    symax = r8vide()
!
    call getvtx(' ', 'TYPE_RESU', scal=typtab, nbret=n1)
!
    call getvtx(' ', 'TYPE_RESU_MECA', scal=typmec, nbret=n1)
!
    call getvtx(' ', 'AXIS', scal=symm, nbret=n1)
    lsymm = .false.
    if (symm .eq. 'OUI') lsymm = .true.
    if (lsymm .and. typmec .ne. 'EVOLUTION') then
        call utmess('F', 'POSTRCCM_59')
    end if
!
!     ------------------------------------------------------------------
!
!     ------------------- TYPE_RESU_MECA = EVOLUTION -------------------
!
!     ------------------------------------------------------------------
!
    if (typmec .eq. 'EVOLUTION') then
!
        call getvtx(' ', 'OPTION', nbval=0, nbret=n1)
        nbopt = -n1
        call getvtx(' ', 'OPTION', nbval=nbopt, vect=kopt, nbret=n1)
!
        call getvid(' ', 'MATER', scal=nommat, nbret=n1)
        call getvr8(' ', 'SY_MAX', scal=symax, nbret=n1)
!
        call rccome(nommat, 'RCCM', icodre)
        if (icodre .eq. 1) then
            call utmess('F', 'POSTRCCM_7', sk='RCCM')
        end if
!
        call rcevol(typtab, nommat, symax, nbopt, kopt, lsymm)
!
!     ------------------------------------------------------------------
!
!     ------------------ TYPE_RESU_MECA = B3600 ------------------
!
!     ------------------------------------------------------------------
!
    else if (typmec .eq. 'B3600') then
!
        call getvtx(' ', 'OPTION', scal=option, nbret=n1)
        if (option .eq. 'FATIGUE') then
            call rc3600()
        elseif (option .eq. 'MOMENT_EQUIVALENT') then
            call rc3600_momeq()
        else
            ASSERT(.false.)
        end if
!
!     ------------------------------------------------------------------
!
!     ---------------TYPE_RESU_MECA = ZE200a, ZE200b, B3200 ----------
!
!     ------------------------------------------------------------------
!
    else
!
        call rc3200()
!
    end if
!
!     ------------------------------------------------------------------
!
    call titre()
!
!
end subroutine
