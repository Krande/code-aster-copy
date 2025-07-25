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
subroutine tbimfi(nparfi, table, newtab, iret)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbextb.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nparfi, iret
    character(len=19) :: table, newtab
!     OPERATEUR  IMPR_TABLE , TRAITEMENT DU MOT CLE FACTEUR "FILTRE"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: jtitr, ititr, ii, ir, ic, ik, ioc, lonmax, lonma1
    integer(kind=8) :: l, l1, l2
    integer(kind=8) :: l3, l4, irt
    character(len=80) :: montit
    character(len=8), pointer :: critere(:) => null()
    character(len=8), pointer :: crit_para(:) => null()
    character(len=24), pointer :: noms_para(:) => null()
    real(kind=8), pointer :: precision(:) => null()
    complex(kind=8), pointer :: vale_c(:) => null()
    integer(kind=8), pointer :: vale_i(:) => null()
    character(len=80), pointer :: vale_k(:) => null()
    real(kind=8), pointer :: vale_r(:) => null()
    character(len=80), pointer :: titr(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    call jeexin(table//'.TITR', irt)
    if (irt .ne. 0) then
        call jeveuo(table//'.TITR', 'L', vk80=titr)
        call jelira(table//'.TITR', 'LONMAX', lonma1)
        lonmax = lonma1+nparfi
        call wkvect(newtab//'.TITR', 'V V K80', lonmax, jtitr)
        do ititr = 1, lonma1
            zk80(jtitr+ititr-1) = titr(ititr)
        end do
    else
        lonma1 = 0
        lonmax = lonma1+nparfi
        call wkvect(newtab//'.TITR', 'V V K80', lonmax, jtitr)
    end if
!
    AS_ALLOCATE(vk24=noms_para, size=nparfi)
    AS_ALLOCATE(vk8=crit_para, size=nparfi)
    AS_ALLOCATE(vi=vale_i, size=nparfi)
    AS_ALLOCATE(vr=vale_r, size=nparfi)
    AS_ALLOCATE(vc=vale_c, size=nparfi)
    AS_ALLOCATE(vk80=vale_k, size=nparfi)
    AS_ALLOCATE(vr=precision, size=nparfi)
    AS_ALLOCATE(vk8=critere, size=nparfi)
!
    ii = -1
    ir = -1
    ic = -1
    ik = -1
!
    do ioc = 1, nparfi
        call getvtx('FILTRE', 'NOM_PARA', iocc=ioc, scal=noms_para(ioc), nbret=l)
        call getvtx('FILTRE', 'CRIT_COMP', iocc=ioc, scal=crit_para(ioc), nbret=l)
        montit = ' '
        call getvis('FILTRE', 'VALE_I', iocc=ioc, nbval=0, nbret=l1)
        call getvr8('FILTRE', 'VALE', iocc=ioc, nbval=0, nbret=l2)
        call getvc8('FILTRE', 'VALE_C', iocc=ioc, nbval=0, nbret=l3)
        call getvtx('FILTRE', 'VALE_K', iocc=ioc, nbval=0, nbret=l4)
        if (l1 .ne. 0) then
            ii = ii+1
            call getvis('FILTRE', 'VALE_I', iocc=ioc, scal=vale_i(ii+1), nbret=l)
            write (montit, 1010) noms_para(ioc), crit_para(ioc), &
                vale_i(ii+1)
        end if
        if (l2 .ne. 0) then
            ir = ir+1
            call getvr8('FILTRE', 'VALE', iocc=ioc, scal=vale_r(ir+1), nbret=l)
            call getvr8('FILTRE', 'PRECISION', iocc=ioc, scal=precision(ir+1), nbret=l)
            call getvtx('FILTRE', 'CRITERE', iocc=ioc, scal=critere(ir+1), nbret=l)
            write (montit, 1020) noms_para(ioc), crit_para(ioc), &
                vale_r(ir+1)
        end if
        if (l3 .ne. 0) then
            ic = ic+1
            call getvc8('FILTRE', 'VALE_C', iocc=ioc, scal=vale_c(ic+1), nbret=l)
            write (montit, 1030) noms_para(ioc), crit_para(ioc), &
                vale_c(ic+1)
        end if
        if (l4 .ne. 0) then
            ik = ik+1
            call getvtx('FILTRE', 'VALE_K', iocc=ioc, scal=vale_k(ik+1), nbret=l)
            write (montit, 1040) noms_para(ioc), crit_para(ioc), &
                vale_k(ik+1)
        end if
        zk80(jtitr+lonma1+ioc-1) = montit
    end do
!
    call tbextb(table, 'V', newtab, nparfi, noms_para, &
                crit_para, vale_i, vale_r, vale_c, vale_k, &
                precision, critere, iret)
!
    AS_DEALLOCATE(vk24=noms_para)
    AS_DEALLOCATE(vk8=crit_para)
    AS_DEALLOCATE(vi=vale_i)
    AS_DEALLOCATE(vr=vale_r)
    AS_DEALLOCATE(vc=vale_c)
    AS_DEALLOCATE(vk80=vale_k)
    AS_DEALLOCATE(vr=precision)
    AS_DEALLOCATE(vk8=critere)
!
1010 format('FILTRE -> NOM_PARA: ', a16, ' CRIT_COMP: ', a4, ' VALE: ', i12)
1020 format('FILTRE -> NOM_PARA: ', a16, ' CRIT_COMP: ', a4,&
   &                                ' VALE: ', 1pe12.5)
1030 format('FILTRE -> NOM_PARA: ', a16, ' CRIT_COMP: ', a4,&
   &                                ' VALE: ', 1pe12.5, 1x, 1pe12.5)
1040 format('FILTRE -> NOM_PARA: ', a16, ' CRIT_COMP: ', a4, ' VALE: ', a8)
!
    call jedema()
end subroutine
