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

subroutine dtmforc_rede(nl_ind, sd_dtm_, sd_nl_, buffdtm, buffnl, &
                        depl, fext)

    implicit none
!
!
! dtmforc_rede : Calculates the force/displacement localized force at
!                the current step (t)
!
!       nl_ind           : nonlinearity index (for sd_nl access)
!       sd_dtm_, buffdtm : dtm data structure and its buffer
!       sd_nl_ , buffnl  : nl  data structure and its buffer
!       depl             : structural modal displacement at instant t
!       fext             : projected total non-linear force
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dtmget.h"
#include "asterfort/fointe.h"
#include "asterfort/nlget.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
!   -0.1- Input/output arguments
    integer(kind=8), intent(in)  :: nl_ind
    character(len=*), intent(in)  :: sd_dtm_
    character(len=*), intent(in)  :: sd_nl_
    integer(kind=8), pointer  :: buffdtm(:)
    integer(kind=8), pointer  :: buffnl(:)
    real(kind=8), pointer  :: depl(:)
    real(kind=8), pointer :: fext(:)
!
!   -0.2- Local variables
    character(len=8)  :: sd_dtm, sd_nl, fonc, comp
    integer(kind=8)           :: im, ier, icomp, nbmode, saredi
    integer(kind=8)           :: start, finish
    real(kind=8)      :: seuil, force
!
    integer(kind=8), pointer :: vindx(:) => null()
    real(kind=8), pointer :: dplred(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    real(kind=8), pointer :: fext0(:) => null()
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
!
    call nlget(sd_nl, _INTERNAL_VARS, vr=vint, buffer=buffnl)
    call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx, buffer=buffnl)
    start = vindx(nl_ind)
!
    call nlget(sd_nl, _FX_FONCT, iocc=nl_ind, kscal=fonc, buffer=buffnl)
    call nlget(sd_nl, _CMP_NAME, iocc=nl_ind, kscal=comp, buffer=buffnl)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode, buffer=buffdtm)
    call nlget(sd_nl, _MODAL_DEPL_NO1, iocc=nl_ind, vr=dplred, buffer=buffnl)
!
    if (comp(1:2) .eq. 'DX') icomp = 1
    if (comp(1:2) .eq. 'DY') icomp = 2
    if (comp(1:2) .eq. 'DZ') icomp = 3
    if (comp(1:3) .eq. 'DRX') icomp = 4
    if (comp(1:3) .eq. 'DRY') icomp = 5
    if (comp(1:3) .eq. 'DRZ') icomp = 6
!
!seuil : deplacement suivant icomp en repere physique
!        (equivalent tophys)
    seuil = 0.d0
    do im = 1, nbmode
        seuil = seuil+dplred((im-1)*6+icomp)*depl(im)
    end do
!
    saredi = 1
!
    call fointe('F ', fonc, 1, [comp], [seuil], &
                force, ier)

    if (abs(force) .le. r8prem()) saredi = 0
!
    AS_ALLOCATE(vr=fext0, size=nbmode)
    fext0(:) = 0.d0

!fext0 : force en repere generalise (equivalent togene)
    do im = 1, nbmode
        fext0(im) = fext0(im)+dplred((im-1)*6+icomp)*force
    end do

    do im = 1, nbmode
        fext(im) = fext(im)+fext0(im)
    end do
!
    AS_DEALLOCATE(vr=fext0)

! --------------------------------------------------------------------------------------------------
!   --- Internal variables, storage
!
    finish = vindx(nl_ind+1)
    ASSERT((finish-start) .eq. NBVARINT_FXRL)

    vint(start) = seuil
    vint(start+1) = force
    vint(start+2) = 1.d0*saredi

end subroutine
