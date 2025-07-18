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
!
subroutine op0042()
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/gettco.h"
#include "asterfort/ccvrpu.h"
#include "asterfort/cresol.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecalr.h"
#include "asterfort/medom1.h"
#include "asterfort/onerrf.h"
#include "asterfort/rsutnu.h"
#include "asterfort/thcalr.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
!
! Command: CALC_ERREUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6) :: nompro
    parameter(nompro='OP0042')
    integer(kind=8) :: ifm, niv, n0, nuord, nchar, ibid, jordr, np, nc
    integer(kind=8) :: nbordr, iret
    real(kind=8) :: prec
    character(len=8) :: resuc1, resuco, modele, cara, crit
    character(len=16) :: nomcmd, tysd, pheno, concep, k16bid, compex
    character(len=19) :: knum, kcha, solveu
    character(len=24) :: mate, mateco
    aster_logical :: newcal
    mpi_int :: msize
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    kcha = '&&'//nompro//'.CHARGES   '
    knum = '&&'//nompro//'.NUME_ORDRE'
!
    call infmaj()
    call infniv(ifm, niv)
!
! --- ON STOCKE LE COMPORTEMENT EN CAS D'ERREUR AVANT MNL : COMPEX
! --- PUIS ON PASSE DANS LE MODE "VALIDATION DU CONCEPT EN CAS D'ERREUR"
    call onerrf(' ', compex, ibid)
    call onerrf('EXCEPTION+VALID', k16bid, ibid)
!
! --- ON INTERDIT L'UTILISATION DE CALC_ERREUR EN PARALLELE
!
    call asmpi_info(size=msize)
    if (msize > 1) then
        call utmess('F', 'CALCULEL3_5', si=to_aster_int(msize))
    end if
!
    call getres(resuc1, concep, nomcmd)
    call getvid(' ', 'RESULTAT', scal=resuco, nbret=n0)
!
    newcal = .false.
    call jeexin(resuc1//'           .DESC', iret)
    if (iret .eq. 0) newcal = .true.
    call gettco(resuco, tysd)
!
    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call rsutnu(resuco, ' ', 0, knum, nbordr, &
                prec, crit, iret)
    if (iret .eq. 10) then
        call utmess('A', 'CALCULEL4_8', sk=resuco)
        goto 999
    end if
    if (iret .ne. 0) then
        call utmess('A', 'ALGORITH3_41')
        goto 999
    end if
!
!     -- ON VEUT INTERDIRE LA REENTRANCE DE LA COMMANDE SI
!        ON UTILISE L'UN DES MOTS CLES : MODELE, CARA_ELEM,
!        CHAM_MATER, EXCIT, GROUP_MA OU MAILLE
!     --------------------------------------------------------
    if (resuco .eq. resuc1) then
        call ccvrpu(resuco, knum, nbordr)
    end if
!
    call jeveuo(knum, 'L', jordr)
    nuord = zi(jordr)
!
!     -- CREATION DU SOLVEUR :
    solveu = '&&OP0042.SOLVEUR'
    call cresol(solveu)
!
    call medom1(modele, mate, mateco, cara, kcha, nchar, &
                resuco, nuord)
    call dismoi('PHENOMENE', modele, 'MODELE', repk=pheno)
!
!     --- TRAITEMENT DU PHENOMENE MECANIQUE ---
    if (pheno(1:4) .eq. 'MECA') then
!
        call mecalr(newcal, tysd, knum, kcha, resuco, &
                    resuc1, nbordr, modele, mate, cara, &
                    nchar)
!
!     --- TRAITEMENT DES PHENOMENES THERMIQUES ET ACOUSTIQUES ---
    else if (pheno(1:4) .eq. 'THER') then
!
        call thcalr(newcal, tysd, knum, kcha, resuco, &
                    resuc1, nbordr, modele, mate, cara, &
                    nchar)
!
    end if
!
999 continue
!
!
! --- ON REMET LE MECANISME D'EXCEPTION A SA VALEUR INITIALE
    call onerrf(compex, k16bid, ibid)

!
    call jedema()
end subroutine
