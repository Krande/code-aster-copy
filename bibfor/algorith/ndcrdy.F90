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
subroutine ndcrdy(result, sddyna)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: result
    character(len=19), intent(out) :: sddyna
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Dynamic management
!
! Create datastructure for dynamic
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of result datastructure (EVOL_NOLI)
! Out sddyna           : name of datastructure for dynamic parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=24) :: tsch, psch, losd, nosd, tfor, cfsc
    integer(kind=8) :: jtsch, jpsch, jlosd, jnosd, jtfor, jcfsc
    character(len=24) :: tcha, ncha, veol, vaol
    integer(kind=8) :: jtcha, jncha, jveol, jvaol
    character(len=24) :: vecent, vecabs
    integer(kind=8) :: jvecen, jvecab
    character(len=16) :: k16bid, nomcmd
    character(len=8) :: k8bid
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)

! --- OPERATEUR APPELANT (STATIQUE OU DYNAMIQUE)
    call getres(k8bid, k16bid, nomcmd)

! - Name of datastructure
    sddyna = result(1:8)//'.SDDYNA'

! - Create objects
    tsch = sddyna(1:15)//'.TYPE_SCH'
    call wkvect(tsch, 'V V K16', 9, jtsch)
    if (nomcmd(1:4) .eq. 'DYNA') then
        zk16(jtsch+1-1) = 'DYNAMIQUE'
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE15_1')
        end if
        psch = sddyna(1:15)//'.PARA_SCH'
        losd = sddyna(1:15)//'.INFO_SD'
        nosd = sddyna(1:15)//'.NOM_SD'
        tfor = sddyna(1:15)//'.TYPE_FOR'
        cfsc = sddyna(1:15)//'.COEF_SCH'
        tcha = sddyna(1:15)//'.TYPE_CHA'
        ncha = sddyna(1:15)//'.NBRE_CHA'
        veol = sddyna(1:15)//'.VEEL_OLD'
        vaol = sddyna(1:15)//'.VEAS_OLD'
        vecent = sddyna(1:15)//'.VECENT'
        vecabs = sddyna(1:15)//'.VECABS'
        call wkvect(psch, 'V V R', 7, jpsch)
        call wkvect(losd, 'V V L', 16, jlosd)
        call wkvect(nosd, 'V V K24', 8, jnosd)
        call wkvect(tfor, 'V V I', 2, jtfor)
        call wkvect(cfsc, 'V V R', 25, jcfsc)
        call wkvect(tcha, 'V V K24', 4, jtcha)
        call wkvect(ncha, 'V V I', 5, jncha)
        call wkvect(veol, 'V V K24', 15, jveol)
        call wkvect(vaol, 'V V K24', 15, jvaol)
        call wkvect(vecent, 'V V K24', 3, jvecen)
        call wkvect(vecabs, 'V V K24', 3, jvecab)
    else
        zk16(jtsch+1-1) = 'STATIQUE'
    end if
!
    call jedema()
!
end subroutine
