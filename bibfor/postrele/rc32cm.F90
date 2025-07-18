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
subroutine rc32cm()
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/codent.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE B3200 et ZE200
!     LECTURE DU MOT CLE FACTEUR "CHAR_MECA"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbchar, iocc, nume, n1, jchar, i, n2
    character(len=8) :: knumec
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
    call getfac('CHAR_MECA', nbchar)
!
    if (nbchar .eq. 0) goto 999
!
    call jecrec('&&RC3200.VALE_CHAR', 'V V R', 'NO', 'DISPERSE', 'VARIABLE', &
                nbchar)
!
    do iocc = 1, nbchar, 1
!
        call getvis('CHAR_MECA', 'NUME_CHAR', iocc=iocc, scal=nume, nbret=n1)
        knumec = 'C       '
        call codent(nume, 'D0', knumec(2:8))
!
        call jecroc(jexnom('&&RC3200.VALE_CHAR', knumec))
        call jeecra(jexnom('&&RC3200.VALE_CHAR', knumec), 'LONMAX', 12)
        call jeecra(jexnom('&&RC3200.VALE_CHAR', knumec), 'LONUTI', 12)
        call jeveuo(jexnom('&&RC3200.VALE_CHAR', knumec), 'E', jchar)
!
!-- cas simple ou cas corps/tubulure ?
        do i = 1, 12
            zr(jchar-1+i) = 0.d0
        end do
!
        call getvr8('CHAR_MECA', 'MX', iocc=iocc, nbval=0, nbret=n2)
!
        if (n2 .ne. 0) then
            call getvr8('CHAR_MECA', 'FX', iocc=iocc, scal=zr(jchar-1+1), nbret=n1)
            call getvr8('CHAR_MECA', 'FY', iocc=iocc, scal=zr(jchar-1+2), nbret=n1)
            call getvr8('CHAR_MECA', 'FZ', iocc=iocc, scal=zr(jchar-1+3), nbret=n1)
            call getvr8('CHAR_MECA', 'MX', iocc=iocc, scal=zr(jchar-1+4), nbret=n1)
            call getvr8('CHAR_MECA', 'MY', iocc=iocc, scal=zr(jchar-1+5), nbret=n1)
            call getvr8('CHAR_MECA', 'MZ', iocc=iocc, scal=zr(jchar-1+6), nbret=n1)
!
        else
            call getvr8('CHAR_MECA', 'FX_TUBU', iocc=iocc, scal=zr(jchar-1+1), nbret=n1)
            call getvr8('CHAR_MECA', 'FY_TUBU', iocc=iocc, scal=zr(jchar-1+2), nbret=n1)
            call getvr8('CHAR_MECA', 'FZ_TUBU', iocc=iocc, scal=zr(jchar-1+3), nbret=n1)
            call getvr8('CHAR_MECA', 'MX_TUBU', iocc=iocc, scal=zr(jchar-1+4), nbret=n1)
            call getvr8('CHAR_MECA', 'MY_TUBU', iocc=iocc, scal=zr(jchar-1+5), nbret=n1)
            call getvr8('CHAR_MECA', 'MZ_TUBU', iocc=iocc, scal=zr(jchar-1+6), nbret=n1)
!
            call getvr8('CHAR_MECA', 'FX_CORP', iocc=iocc, scal=zr(jchar-1+7), nbret=n1)
            call getvr8('CHAR_MECA', 'FY_CORP', iocc=iocc, scal=zr(jchar-1+8), nbret=n1)
            call getvr8('CHAR_MECA', 'FZ_CORP', iocc=iocc, scal=zr(jchar-1+9), nbret=n1)
            call getvr8('CHAR_MECA', 'MX_CORP', iocc=iocc, scal=zr(jchar-1+10), nbret=n1)
            call getvr8('CHAR_MECA', 'MY_CORP', iocc=iocc, scal=zr(jchar-1+11), nbret=n1)
            call getvr8('CHAR_MECA', 'MZ_CORP', iocc=iocc, scal=zr(jchar-1+12), nbret=n1)
        end if
!
    end do
!
999 continue
    call jedema()
end subroutine
