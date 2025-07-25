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

subroutine op0052()
! person_in_charge: nicolas.sellenet at edf.fr
! ----------------------------------------------------------------------
!  COMMANDE CALC_CHAMP
! ----------------------------------------------------------------------
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/ccbcop.h"
#include "asterfort/ccchut.h"
#include "asterfort/cclopu.h"
#include "asterfort/ccvrpu.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/onerrf.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
    character(len=6) :: nompro
    parameter(nompro='OP0052')
!
    integer(kind=8) :: ifm, niv, ibid, n0, iret, np, nc
    integer(kind=8) :: nbordr, nbropt
!
    real(kind=8) :: prec
!
    character(len=8) :: resuc1, resuco, crit
    character(len=16) :: compex, k16bid, concep, nomcmd
    character(len=19) :: lisord, lisopt
!
    call jemarq()
!
    lisopt = '&&'//nompro//'.LIS_OPTION'
    lisord = '&&'//nompro//'.NUME_ORDRE'
!
    call infmaj()
    call infniv(ifm, niv)
!
!     ON STOCKE LE COMPORTEMENT EN CAS D'ERREUR AVANT MNL : COMPEX
!     PUIS ON PASSE DANS LE MODE "VALIDATION DU CONCEPT EN CAS D'ERREUR"
    call onerrf(' ', compex, ibid)
    call onerrf('EXCEPTION+VALID', k16bid, ibid)
!
!     RECUPERATION DES NOMS DES SD RESULTAT
    call getres(resuc1, concep, nomcmd)
    call getvid(' ', 'RESULTAT', scal=resuco, nbret=n0)
!
    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call rsutnu(resuco, ' ', 0, lisord, nbordr, &
                prec, crit, iret)
    if (iret .eq. 10) then
        call utmess('A', 'CALCULEL4_8', sk=resuco)
        goto 9999
    end if
    if (iret .ne. 0) then
        call utmess('A', 'ALGORITH3_41')
        goto 9999
    end if
!
!     ON VEUT INTERDIRE LA REENTRANCE DE LA COMMANDE SI
!     ON UTILISE L'UN DES MOTS CLES : MODELE, CARAEL_ELEM,
!     CHAM_CHMATER OU EXCIT
    if (resuco .eq. resuc1) then
        call ccvrpu(resuco, lisord, nbordr)
    end if
!
!     FABRICATION DE LA LISTE DES OPTIONS
    call cclopu(resuco, resuc1, lisord, nbordr, lisopt, &
                nbropt)
!
!     APPEL A LA ROUTINE PREPARANT L'APPEL A CALCOP
    call ccbcop(resuco, resuc1, lisord, nbordr, lisopt, &
                nbropt)
!
    call jedetr(lisopt)
!
    call ccchut(resuco, resuc1, lisord, nbordr)
!
9999 continue
!     ON REMET LE MECANISME D'EXCEPTION A SA VALEUR INITIALE
    call onerrf(compex, k16bid, ibid)
!
    call refdcp(resuco, resuc1)
!
    call jedema()
!
end subroutine
