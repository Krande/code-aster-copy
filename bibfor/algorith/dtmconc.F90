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

subroutine dtmconc(sd_dtm_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmconc : Append several dyna_gene data structures saved by an adaptative time-step
!           integration algorithm into a single result.
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dtmget.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/nlget.h"
#include "asterfort/wkvect.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
!
!   -0.2- Local variables
    integer(kind=8)           :: iarch_sd, i_bloc, nbmode, nbnoli, nbvint, l_disc, l_ordr, i_disc
    integer(kind=8)           :: ordr, ordr_prec
    real(kind=8)      :: disc, disc_prec
    character(len=7)  :: intk7
    character(len=8)  :: sd_dtm
    character(len=8)  :: nomres, sd_nl
    character(len=16) :: nomres16
    integer(kind=8), pointer  :: isto(:) => null()
    integer(kind=8), pointer :: vindx(:) => null()
    real(kind=8), pointer  :: v_bloc(:) => null()
    real(kind=8), pointer  :: v_disc(:) => null()
    integer(kind=8), pointer  :: v_blo2(:) => null()
    integer(kind=8), pointer  :: v_ordr(:) => null()
!
!   0 - Initializations
    call jemarq()

    sd_dtm = sd_dtm_

    call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
    call dtmget(sd_dtm, _IARCH_SD, iscal=iarch_sd)
    call dtmget(sd_dtm, _ARCH_STO, vi=isto)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnoli)
    if (nbnoli .gt. 0) then
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)
    end if

    call wkvect(nomres//'           .BLOC', 'G V R', iarch_sd, vr=v_bloc)
    call wkvect(nomres//'           .BLO2', 'G V I', iarch_sd, vi=v_blo2)

    disc_prec = -1
    ordr_prec = -1
    do i_bloc = 1, iarch_sd
        call codent(i_bloc, 'D0', intk7)
        nomres16 = nomres//'.'//intk7
        if (i_bloc .eq. iarch_sd) then
            call juveca(nomres16//'   .ORDR', isto(1))
            call juveca(nomres16//'   .DISC', isto(1))
            call juveca(nomres16//'   .PTEM', isto(1))
            call juveca(nomres16//'   .DEPL', isto(1)*nbmode)
            call juveca(nomres16//'   .VITE', isto(1)*nbmode)
            call juveca(nomres16//'   .ACCE', isto(1)*nbmode)
!           Nonlinearities
            if (nbnoli .gt. 0) then
                call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx)
                nbvint = vindx(nbnoli+1)-1
                call juveca(nomres16//'.NL.VINT', isto(1)*nbvint)
            end if
        end if
        call jeveuo(nomres16//'   .DISC', 'L', vr=v_disc)
        call jelira(nomres16//'   .DISC', 'LONMAX', l_disc)
        call jeveuo(nomres16//'   .ORDR', 'L', vi=v_ordr)
        call jelira(nomres16//'   .ORDR', 'LONMAX', l_ordr)
        ASSERT(l_disc .eq. l_ordr)
        do i_disc = 1, l_disc
            disc = v_disc(i_disc)
            ordr = v_ordr(i_disc)
            ASSERT(disc_prec .le. disc .and. ordr_prec .le. ordr)
            disc_prec = disc
            ordr_prec = ordr
        end do
        v_bloc(i_bloc) = disc_prec
        v_blo2(i_bloc) = ordr_prec
        call jelibe(nomres16//'   .DISC')
        call jelibe(nomres16//'   .ORDR')
    end do
!
    call jedema()
end subroutine
