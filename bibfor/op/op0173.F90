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
subroutine op0173()
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/putvir.h"
#include "asterc/putvrr.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/sdmpic.h"
#include "asterfort/tbcopi.h"
#include "asterfort/tbimfi.h"
#include "asterfort/tbliva.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
! ----------------------------------------------------------------------
!
!
!     COMMANDE:  EXTR_TABLE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, n1, iret, nparfi, vali
    real(kind=8) :: r8b, valr
    complex(kind=8) :: cbid, valc
    character(len=8) :: k8b, nomres, ctype, table
    character(len=8) :: nomgd, tsca
    character(len=16) :: nomcmd, concep, typesd
    character(len=19) :: newtab, newta1
    character(len=24) :: para
    character(len=80) :: valk
!     ------------------------------------------------------------------
!
    call jemarq()
    cbid = (0.d0, 0.d0)
    r8b = 0.d0
    ibid = 0
!
    call getres(nomres, concep, nomcmd)
!
    call getvid(' ', 'TABLE', scal=table, nbret=n1)
    newtab = table
!
    call getvtx(' ', 'NOM_PARA', scal=para, nbret=n1)
!
    call getvtx(' ', 'TYPE_RESU', scal=typesd, nbret=n1)
!
    call getfac('FILTRE', nparfi)
    if (nparfi .ne. 0) then
        newta1 = '&&OP0173.FILTRE '
        call tbimfi(nparfi, newtab, newta1, iret)
        if (iret .ne. 0) then
            call utmess('F', 'UTILITAI7_11')
        end if
        newtab = newta1
    end if
!
    call tbliva(newtab, 0, k8b, [ibid], [r8b], &
                [cbid], k8b, k8b, [r8b], para, &
                ctype, vali, valr, valc, valk, &
                iret)
    if (iret .eq. 0) then
    else if (iret .eq. 1) then
        call utmess('F', 'CALCULEL4_43')
    else if (iret .eq. 2) then
        call utmess('F', 'CALCULEL4_44')
    else if (iret .eq. 3) then
        call utmess('F', 'CALCULEL4_45')
    else
        call utmess('F', 'CALCULEL4_46')
    end if
!
    if (typesd .eq. 'MATR_ASSE_GENE_R') then
!          ------------------------------
        call copisd('MATR_ASSE_GENE', 'G', valk, nomres)
!
    else if (typesd .eq. 'MATR_ELEM_DEPL_R') then
!          ------------------------------
        call copisd('MATR_ELEM', 'G', valk, nomres)
        call sdmpic('MATR_ELEM', nomres)
!
    else if (typesd .eq. 'VECT_ELEM_DEPL_R') then
!          ------------------------------
        call copisd('VECT_ELEM', 'G', valk, nomres)
!
    else if (typesd .eq. 'MATR_ELEM_TEMP_R') then
!          ------------------------------
        call copisd('MATR_ELEM', 'G', valk, nomres)
        call sdmpic('MATR_ELEM', nomres)
!
    else if (typesd .eq. 'VECT_ELEM_TEMP_R') then
!          ------------------------------
        call copisd('VECT_ELEM', 'G', valk, nomres)
!
    elseif (typesd .eq. 'CHAM_GD_SDASTER' .or. typesd .eq. &
            'CHAM_NO_SDASTER' .or. typesd .eq. 'CARTE_SDASTER' .or. typesd &
            .eq. 'CHAM_ELEM') then
!          ----------------------------------------
        call dismoi('NOM_GD', valk, 'CHAMP', repk=nomgd)
        call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
        if (tsca .ne. 'R') then
            call utmess('F', 'CALCULEL4_27')
        end if
        call copisd('CHAMP_GD', 'G', valk, nomres)
!
    else if (typesd .eq. 'MODE_MECA') then
!          ------------------------------
        call copisd('RESULTAT', 'G', valk, nomres)
!
    elseif (typesd .eq. 'FONCTION_SDASTER' .or. typesd .eq. &
            'FONCTION_C' .or. typesd .eq. 'NAPPE_SDASTER') then
!          ------------------------------
        call copisd('FONCTION', 'G', valk, nomres)
!
    else if (typesd .eq. 'TABLE_SDASTER') then
        call tbcopi('G', valk, nomres)
!
    else if (typesd .eq. 'ENTIER') then
        call putvir(vali)
!
    else if (typesd .eq. 'REEL') then
        call putvrr(valr)
!
    else
        call utmess('F', 'CALCULEL4_47', sk=typesd)
    end if
!
    if (typesd .eq. 'REEL' .and. typesd .eq. 'ENTIER') then
        call titre()
    end if
!
    call jedema()
!
end subroutine
