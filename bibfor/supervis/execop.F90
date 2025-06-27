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

subroutine execop(num)
    use superv_module, only: superv_before, superv_after
    implicit none
!     EXECUTION DE LA COMMANDE
!     ------------------------------------------------------------------
!     COMMON POUR LE NIVEAU D'"INFO"
#include "asterc/etausr.h"
#include "asterc/gcecdu.h"
#include "asterc/uttrst.h"
#include "asterf.h"
#include "asterfort/assert.h"
#include "asterfort/ex0000.h"
#include "asterfort/iunifi.h"
#include "asterfort/jevema.h"
#include "asterfort/sigusr.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpg.h"
    integer(kind=8), intent(in), optional :: num
!
    integer(kind=8) :: nivuti, nivpgm, unite
    common/inf001/nivuti, nivpgm, unite
!
    integer(kind=8) :: nuoper, imaav, imaap
    real(kind=8) :: tpres
!     ------------------------------------------------------------------
!
    if (present(num)) then
        nuoper = num
    else
        call gcecdu(nuoper)
    end if
!
    if (nuoper .le. 0) then
        ! only DEBUT in this case
        ASSERT(nuoper .eq. 0)
        call ex0000(nuoper)
    else
!
!     -- ON NOTE LA MARQUE AVANT D'APPELER LA PROCHAINE COMMANDE :
        call jevema(imaav)
!
        call superv_before()
!
!     -- ON INITIALISATION DES COMPTEURS DE TEMPS :
        call uttcpg('INIT', ' ')
!
!     -- ON MET A JOUR LE COMMON INF001 :
        nivuti = 1
        nivpgm = 1
        unite = iunifi('MESSAGE')
!
        if (nuoper .lt. 200) then
            call ex0000(nuoper)
        else if (nuoper .eq. 8888) then
!       special operator number: does nothing, just to pass after/before steps
        else
            call utmess('E', 'SUPERVIS_61', si=nuoper)
        end if
!
! --- VERIFICATION SI INTERRUPTION DEMANDEE PAR SIGNAL USR1
!
        if (etausr() .eq. 1) then
            call sigusr()
        end if
!
        call uttrst(tpres)
        if (tpres .lt. 0.d0) then
            call utmess('Z', 'SUPERVIS_63', sr=-tpres, num_except=ASTER_TIMELIMIT_ERROR)
        end if
!
!     -- CONTROLE DE L'APPARIEMMENT DES JEMARQ/JEDEMA
        call jevema(imaap)
        if (imaav .ne. imaap) then
            call utmess('F', 'SUPERVIS_3', sk='JEMARQ/JEDEMA')
        end if
!
!     -- ON IMPRIME LES COMPTEURS DE TEMPS :
!        (IL FAUT LE FAIRE AVANT LA DESTRUCTION DES OBJETS VOLATILES)
        call uttcpg('IMPR', 'CUMU')
!
        call superv_after()
    end if
!
end subroutine
