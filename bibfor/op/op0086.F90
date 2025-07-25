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

subroutine op0086()
    implicit none
!     COMMANDE:  MACR_ELEM_STAT
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sschge.h"
#include "asterfort/ssdege.h"
#include "asterfort/ssmage.h"
#include "asterfort/ssrige.h"
#include "asterfort/utmess.h"
    character(len=8) :: nomu, macrelem_reuse
    character(len=16) :: kbi1, kbi2
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iarefm, iret, nocc
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    call getres(nomu, kbi1, kbi2)

    call getvid(' ', 'MACR_ELEM', scal=macrelem_reuse, nbret=nocc)
    if (nocc .ne. 0) then
        if (nomu .ne. macrelem_reuse) then
            call utmess('F', 'SUPERVIS2_79', sk='MACR_ELEM')
        end if
    end if
!
!
!     --TRAITEMENT DES MOTS CLEFS 'DEFINITION' ET 'EXTERIEUR'
!     --------------------------------------------------------
    call getfac('DEFINITION', nocc)
    if (nocc .eq. 1) then
        call jeexin(nomu//'.REFM', iret)
        if (iret .gt. 0) then
            call utmess('F', 'SOUSTRUC_9', sk=nomu)
        else
            call ssdege(nomu)
        end if
    end if
!
!
!     --TRAITEMENT DU MOT CLEF 'RIGI_MECA'
!     --------------------------------------
    call getfac('RIGI_MECA', nocc)
    if (nocc .eq. 1) then
        call jeexin(nomu//'.REFM', iret)
        if (iret .eq. 0) then
            call utmess('F', 'SOUSTRUC_10')
        end if
        call jeveuo(nomu//'.REFM', 'L', iarefm)
        if (zk8(iarefm-1+6) .eq. 'OUI_RIGI') then
            call utmess('F', 'SOUSTRUC_11', sk=nomu)
        else
            call ssrige(nomu)
        end if
    end if
!
!
!     --TRAITEMENT DU MOT CLEF 'MASS_MECA':
!     --------------------------------------
    call getfac('MASS_MECA', nocc)
    if (nocc .eq. 1) then
        call jeexin(nomu//'.REFM', iret)
        if (iret .eq. 0) then
            call utmess('F', 'SOUSTRUC_12')
        end if
        call jeveuo(nomu//'.REFM', 'L', iarefm)
        if (zk8(iarefm-1+6) .ne. 'OUI_RIGI') then
            call utmess('F', 'SOUSTRUC_12')
        end if
        if (zk8(iarefm-1+7) .eq. 'OUI_MASS') then
            call utmess('F', 'SOUSTRUC_13')
        else
            call ssmage(nomu, 'MASS_MECA')
        end if
    end if
!
!
!     --TRAITEMENT DU MOT CLEF 'AMOR_MECA':
!     --------------------------------------
    call getfac('AMOR_MECA', nocc)
    if (nocc .eq. 1) then
        call jeexin(nomu//'.REFM', iret)
        if (iret .eq. 0) then
            call utmess('F', 'SOUSTRUC2_4')
        end if
        call jeveuo(nomu//'.REFM', 'L', iarefm)
        if (zk8(iarefm-1+8) .eq. 'OUI_AMOR') then
            call utmess('F', 'SOUSTRUC2_5', sk=nomu)
        else
            call ssmage(nomu, 'AMOR_MECA')
        end if
    end if
!
!
!     --TRAITEMENT DU MOT CLEF 'CAS_CHARGE'
!     --------------------------------------
    call getfac('CAS_CHARGE', nocc)
    if (nocc .gt. 0) then
        call jeexin(nomu//'.REFM', iret)
        if (iret .eq. 0) then
            call utmess('F', 'SOUSTRUC_14')
        end if
        call jeveuo(nomu//'.REFM', 'L', iarefm)
        if (zk8(iarefm-1+6) .ne. 'OUI_RIGI') then
            call utmess('F', 'SOUSTRUC_14')
        end if
        call sschge(nomu)
    end if
!
!
    call jedema()
end subroutine
