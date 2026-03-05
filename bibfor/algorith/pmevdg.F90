! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine pmevdg(sddisc, tablType, tablIncr, iEvent, iEventActi)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbacce.h"
#include "asterfort/tbliva.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "event_def.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: tablType
    character(len=24), intent(in) :: tablIncr
    integer(kind=8), intent(in) :: iEvent
    integer(kind=8), intent(out) :: iEventActi
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - EVENEMENTS)
!
! GESTION DE L'EVENEMENT DELTA_GRANDEUR
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, ier
    integer(kind=8) :: etat_loca
    integer(kind=8), pointer:: loca(:) => null()
    real(kind=8) :: valeRefe, vale, r8bid
    integer(kind=8) :: ibid
    character(len=8) :: k8bid, crit
    complex(kind=8) :: c16bid
    character(len=16) :: cmpName
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<SIMUPOINTMAT> ... DELTA_GRANDEUR'
    end if

! - INITIALISATIONS
    iEventActi = 0
    r8bid = 0.d0

! - PARAMETRES
    call utdidt('L', sddisc, 'ECHE', 'NOM_CMP', index_=iEvent, valk_=cmpName)
    call utdidt('L', sddisc, 'ECHE', 'VALE_REF', index_=iEvent, valr_=valeRefe)
    call utdidt('L', sddisc, 'ECHE', 'CRIT_COMP', index_=iEvent, valk_=crit)
    call jeveuo(sddisc//'.ELOC', 'L', vi=loca)
    etat_loca = loca(SIZE_LELOCA*(iEvent-1)+1)
    if (etat_loca .ne. LOCA_TOUT) then
        call utmess('F', 'COMPOR2_9')
    end if

! - Get value
    if (tablType .eq. 0) then
        call tbacce(tablIncr, 1, cmpName, 'L', ibid, vale, c16bid, k8bid)

    else if (tablType .eq. 1) then
        call tbliva(tablIncr, 1, 'CMP', [ibid], [r8bid], &
                    [c16bid], cmpName, 'EGAL', [0.d0], 'VALEUR', &
                    k8bid, ibid, vale, c16bid, k8bid, &
                    ier)
        if (ier .ne. 0) then
            vale = 0.d0
        end if

    else
        ASSERT(ASTER_FALSE)
    end if
!
    ASSERT(crit .eq. 'GT')
!
    if (vale .gt. valeRefe) then
        iEventActi = iEvent
    end if
!
    call jedema()
end subroutine
