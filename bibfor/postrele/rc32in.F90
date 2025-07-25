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
subroutine rc32in()
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE B3200 et ZE200
!     RECUPERATION
!          DE  C1, C2, K1, K2, C3, K3      SOUS INDI_SIGM
!          DE  DIAM, INERTIE, EP           SOUS CARAC_TUYAU
!          DE  KT_SN, KT_SP                SOUS FACT_SIGM
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ndim, jvalin, n1, ktsn, ktsp, n2, n2a, n2b, n2c, n2d
    real(kind=8) :: bid
    integer(kind=8) :: n2e, n2f, n2g, n2h, i, j, n1a, n1b, n1c
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
    ndim = 19
    call wkvect('&&RC3200.INDI', 'V V R', ndim, jvalin)
!
! --- indices de contrainte et carac. géométriques de la tuyauterie
! ----- cas du ze200
    call getvr8('INDI_SIGM', 'K1', nbval=0, iocc=1, nbret=n1)
    if (n1 .ne. 0) then
        call getvr8('INDI_SIGM', 'K1', scal=zr(jvalin), iocc=1, nbret=n1)
        call getvr8('INDI_SIGM', 'C1', scal=zr(jvalin+1), iocc=1, nbret=n1)
        call getvr8('INDI_SIGM', 'K3', scal=zr(jvalin+4), iocc=1, nbret=n1)
        call getvr8('INDI_SIGM', 'C3', scal=zr(jvalin+5), iocc=1, nbret=n1)
!
        call getvr8('TUYAU', 'R', scal=zr(jvalin+6), iocc=1, nbret=n1)
        call getvr8('TUYAU', 'EP', scal=zr(jvalin+7), iocc=1, nbret=n1)
        do i = 1, 8
            zr(jvalin+i+10) = 1
        end do
        zr(jvalin+15) = 0
        zr(jvalin+17) = 0
!------ cas corps-tubulure
        call getvr8('CHAR_MECA', 'MX_TUBU', iocc=1, scal=bid, nbret=n2)
        if (n2 .ne. 0) then
            call getvr8('INDI_SIGM', 'K2_TUBU', scal=zr(jvalin+11), iocc=1, nbret=n2a)
            call getvr8('INDI_SIGM', 'C2_TUBU', scal=zr(jvalin+12), iocc=1, nbret=n2b)
            call getvr8('INDI_SIGM', 'K2_CORP', scal=zr(jvalin+13), iocc=1, nbret=n2c)
            call getvr8('INDI_SIGM', 'C2_CORP', scal=zr(jvalin+14), iocc=1, nbret=n2d)
            call getvr8('TUYAU', 'R_TUBU', scal=zr(jvalin+15), iocc=1, nbret=n2e)
            call getvr8('TUYAU', 'I_TUBU', scal=zr(jvalin+16), iocc=1, nbret=n2f)
            call getvr8('TUYAU', 'R_CORP', scal=zr(jvalin+17), iocc=1, nbret=n2g)
            call getvr8('TUYAU', 'I_CORP', scal=zr(jvalin+18), iocc=1, nbret=n2h)
            if (n2a*n2b*n2c*n2d*n2e*n2f*n2g*n2h .eq. 0) then
                call utmess('F', 'POSTRCCM_46')
            end if
            call getvr8('INDI_SIGM', 'K2', scal=zr(jvalin+2), iocc=1, nbret=n1a)
            call getvr8('INDI_SIGM', 'C2', scal=zr(jvalin+3), iocc=1, nbret=n1b)
            call getvr8('TUYAU', 'I', scal=zr(jvalin+8), iocc=1, nbret=n1c)
            if (n1a+n1b+n1c .ne. 0) then
                call utmess('F', 'POSTRCCM_46')
            else
                zr(jvalin+2) = 0
                zr(jvalin+3) = 0
                zr(jvalin+3) = 1
            end if
        else
            call getvr8('INDI_SIGM', 'K2', scal=zr(jvalin+2), iocc=1, nbret=n1a)
            call getvr8('INDI_SIGM', 'C2', scal=zr(jvalin+3), iocc=1, nbret=n1b)
            call getvr8('TUYAU', 'I', scal=zr(jvalin+8), iocc=1, nbret=n1c)
            if (n1a*n1b*n1c .eq. 0) call utmess('F', 'POSTRCCM_46')
        end if
! ----- cas du b3200 sans indices de contraintes
    else
        do j = 1, 19
            zr(jvalin+j-1) = 1
        end do
        zr(jvalin+15) = 0
        zr(jvalin+17) = 0
    end if
!
! --- facteur de concentration de contrainte (b3200 uniquement)
    zr(jvalin+9) = 1
    zr(jvalin+10) = 1
!
    call getvr8('FACT_SIGM', 'KT_SN', iocc=1, nbval=0, nbret=ktsn)
    if (ktsn .ne. 0) then
        call getvr8('FACT_SIGM', 'KT_SN', scal=zr(jvalin+9), iocc=1, nbret=ktsn)
    end if
!
    call getvr8('FACT_SIGM', 'KT_SP', iocc=1, nbval=0, nbret=ktsp)
    if (ktsp .ne. 0) then
        call getvr8('FACT_SIGM', 'KT_SP', scal=zr(jvalin+10), iocc=1, nbret=ktsp)
    end if
!
    call jedema()
end subroutine
