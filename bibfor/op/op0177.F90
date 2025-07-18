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

subroutine op0177()
!
!     COMMANDE:  TEST_TABLE
!
! ----------------------------------------------------------------------
    implicit none
!
! 0.1. ==> ARGUMENTS
!
!
! 0.2. ==> COMMUNS
! 0.3. ==> VARIABLES LOCALES
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tbimfi.h"
#include "asterfort/tbliva.h"
#include "asterfort/ulexis.h"
#include "asterfort/ulopen.h"
#include "asterfort/tresu_ordgrd.h"
#include "asterfort/tresu_print_all.h"
#include "asterfort/tresu_tabl.h"
#include "asterfort/tresu_read_refe.h"
#include "asterfort/tresu_str.h"
#include "asterfort/tresu_tole.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=6) :: nompro
    parameter(nompro='OP0177')
!
    integer(kind=8) :: ibid, n1, n2, n3, nk, iret, ific, nparfi
    integer(kind=8) :: vali, irefr, irefi, irefc, irefk, nref
    integer(kind=8) :: n1r, n2r, n3r, nkr, irefrr, irefir, irefcr, irefkr
    integer(kind=8) :: nl1, nl2, nl11, nl22
!
    real(kind=8) :: r8b, valr, epsi, epsir
!
    complex(kind=8) :: cbid, valc
!
    character(len=1) :: typr
    character(len=3) :: ssigne
    character(len=8) :: k8b, crit, ctype, typtes, latabl
    character(len=8) :: motcle
    character(len=16) :: nomfi, tbtxt(2), tbref(2)
    character(len=19) :: newtab, newta1
    character(len=24) :: para
    character(len=24) :: travr, travi, travc, travrr, travir, travcr, travk, travkr
    character(len=80) :: valk
    character(len=200) :: lign1, lign2
    aster_logical :: lref
    aster_logical :: skip
    real(kind=8) :: ordgrd
!     ------------------------------------------------------------------
!
!====
! 1. PREALABLES
!====
!
    call jemarq()
    ibid = 0
    cbid = (0.d0, 0.d0)
    r8b = 0.d0
    call infmaj()
!
    travr = '&&'//nompro//'_TRAVR          '
    travi = '&&'//nompro//'_TRAVI          '
    travc = '&&'//nompro//'_TRAVC          '
    travk = '&&'//nompro//'_TRAVK          '
    travrr = '&&'//nompro//'_TRAVR_R        '
    travir = '&&'//nompro//'_TRAVI_R        '
    travcr = '&&'//nompro//'_TRAVC_R        '
    travkr = '&&'//nompro//'_TRAVK_R        '
    motcle = 'TABLE'
!
    nomfi = ' '
    ific = iunifi('RESULTAT')
    if (.not. ulexis(ific)) then
        call ulopen(ific, ' ', nomfi, 'NEW', 'O')
    end if
    write (ific, 100)
!
    call getvid(' ', 'TABLE', scal=latabl, nbret=n1)
!
    call getfac('FILTRE', nparfi)
!
    call getvtx(' ', 'VALE_ABS', scal=ssigne, nbret=n1)
    call tresu_tole(epsi)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
!
    call getvr8(' ', 'VALE_CALC', nbval=0, nbret=n1)
    call getvis(' ', 'VALE_CALC_I', nbval=0, nbret=n2)
    call getvc8(' ', 'VALE_CALC_C', nbval=0, nbret=n3)
    call getvtx(' ', 'VALE_CALC_K', nbval=0, nbret=nk)
!
    irefr = 0
    irefi = 0
    irefc = 0
    irefk = 0
    skip = .false.
    ordgrd = 1.d0
    if (n1 .ne. 0) then
        nref = -n1
        typr = 'R'
        call jedetr(travr)
        call wkvect(travr, 'V V R', nref, irefr)
        call getvr8(' ', 'VALE_CALC', nbval=nref, vect=zr(irefr))
        call tresu_ordgrd(zr(irefr), skip, ordgrd)
    else if (n2 .ne. 0) then
        nref = -n2
        typr = 'I'
        call jedetr(travi)
        call wkvect(travi, 'V V I', nref, irefi)
        call getvis(' ', 'VALE_CALC_I', nbval=nref, vect=zi(irefi))
    else if (n3 .ne. 0) then
        nref = -n3
        typr = 'C'
        call jedetr(travc)
        call wkvect(travc, 'V V C', nref, irefc)
        call getvc8(' ', 'VALE_CALC_C', nbval=nref, vect=zc(irefc))
    else
        ASSERT(nk .ne. 0)
        nref = -nk
        typr = 'K'
        call jedetr(travk)
        call wkvect(travk, 'V V K80', nref, irefk)
        call getvtx(' ', 'VALE_CALC_K', nbval=nref, vect=zk80(irefk))
    end if
! ----------------------------------------------------------------------
    lref = .false.
    call getvr8(' ', 'PRECISION', scal=epsir, nbret=iret)
    if (iret .ne. 0) then
        lref = .true.
        call getvr8(' ', 'VALE_REFE', nbval=0, nbret=n1r)
        call getvis(' ', 'VALE_REFE_I', nbval=0, nbret=n2r)
        call getvc8(' ', 'VALE_REFE_C', nbval=0, nbret=n3r)
        call getvtx(' ', 'VALE_REFE_K', nbval=0, nbret=nkr)
!
        irefrr = 0
        irefir = 0
        irefcr = 0
        irefkr = 0
        if (n1r .ne. 0) then
            ASSERT((n1r .eq. n1))
            nref = -n1r
            call jedetr(travrr)
            call wkvect(travrr, 'V V R', nref, irefrr)
            call getvr8(' ', 'VALE_REFE', nbval=nref, vect=zr(irefrr))
        else if (n2r .ne. 0) then
            ASSERT((n2r .eq. n2))
            nref = -n2r
            call jedetr(travir)
            call wkvect(travir, 'V V I', nref, irefir)
            call getvis(' ', 'VALE_REFE_I', nbval=nref, vect=zi(irefir))
        else if (n3r .ne. 0) then
            ASSERT((n3r .eq. n3))
            nref = -n3r
            call jedetr(travcr)
            call wkvect(travcr, 'V V C', nref, irefcr)
            call getvc8(' ', 'VALE_REFE_C', nbval=nref, vect=zc(irefcr))
        else
            ASSERT(nkr .ne. 0)
            ASSERT(nkr .eq. nk)
            nref = -nkr
            call jedetr(travkr)
            call wkvect(travkr, 'V V K80', nref, irefkr)
            call getvtx(' ', 'VALE_REFE_K', nbval=nref, vect=zk80(irefkr))
        end if
    end if
    if (skip .and. .not. lref) then
        call utmess('I', 'TEST0_11')
    end if
! ----------------------------------------------------------------------
    call getvtx(' ', 'NOM_PARA', scal=para)
    call getvtx(' ', 'TYPE_TEST', scal=typtes, nbret=n1)
!
    lign1 = ' '
    lign2 = ' '
!
!   Traitement du mot clé "FILTRE" ---
    newtab = latabl
!
    lign1(1:21) = '---- '//motcle
    lign1(22:22) = '.'
    lign2(1:21) = '     '//latabl
    lign2(22:22) = '.'
    nl1 = lxlgut(lign1)
    nl2 = lxlgut(lign2)
    lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_PARA'
    lign2(1:nl2+16) = lign2(1:nl2-1)//' '//para(1:16)
    lign1(nl1+17:nl1+17) = '.'
    lign2(nl2+17:nl2+17) = '.'
!
    if (nparfi .ne. 0) then
        newta1 = '&&'//nompro//'.FILTRE '
        call tbimfi(nparfi, newtab, newta1, iret)
        if (iret .ne. 0) then
            call utmess('F', 'CALCULEL6_7')
        end if
        newtab = newta1
    end if
!   ------------------------------------------------------------------
!
    call tresu_read_refe(' ', 1, tbtxt)
    if (lref) then
        tbref(1) = tbtxt(1)
        tbref(2) = tbtxt(2)
        tbtxt(1) = 'NON_REGRESSION'
    end if
!
    if (n1 .ne. 0) then
!   cas de TYPE_TEST
!
        nl1 = lxlgut(lign1)
        nl2 = lxlgut(lign2)
        lign1(1:nl1+16) = lign1(1:nl1-1)//' TYPE_TEST'
        lign2(1:nl2+16) = lign2(1:nl2-1)//' '//typtes
        lign1(nl1+17:nl1+17) = '.'
        lign2(nl2+17:nl2+17) = '.'
!
        nl1 = lxlgut(lign1)
        nl11 = lxlgut(lign1(1:nl1-1))
        nl2 = lxlgut(lign2)
        nl22 = lxlgut(lign2(1:nl2-1))
        if (nl11 .lt. 80) then
            write (ific, *) lign1(1:nl11)
        else if (nl11 .lt. 160) then
            write (ific, 116) lign1(1:80), lign1(81:nl11)
        else
            write (ific, 120) lign1(1:80), lign1(81:160), lign1(161: &
                                                                nl11)
        end if
        if (nl22 .lt. 80) then
            write (ific, *) lign2(1:nl22)
        else if (nl22 .lt. 160) then
            write (ific, 116) lign2(1:80), lign2(81:nl22)
        else
            write (ific, 120) lign2(1:80), lign2(81:160), lign2(161: &
                                                                nl22)
        end if
!
        call tresu_tabl(newtab, para, typtes, typr, tbtxt, &
                        zi(irefi), zr(irefr), zc(irefc), epsi, crit, &
                        .true._1, ssigne, ignore=skip, compare=ordgrd)
        if (lref) then
            call tresu_tabl(newtab, para, typtes, typr, tbref, &
                            zi(irefir), zr(irefrr), zc(irefcr), epsir, crit, &
                            .false._1, ssigne)
        end if
    else
!
        call tbliva(newtab, 0, k8b, [ibid], [r8b], &
                    [cbid], k8b, k8b, [r8b], para, &
                    ctype, vali, valr, valc, valk, &
                    iret)
!
!
        nl1 = lxlgut(lign1)
        nl11 = lxlgut(lign1(1:nl1-1))
        nl2 = lxlgut(lign2)
        nl22 = lxlgut(lign2(1:nl2-1))
        if (nl11 .lt. 80) then
            write (ific, *) lign1(1:nl11)
        else if (nl11 .lt. 160) then
            write (ific, 116) lign1(1:80), lign1(81:nl11)
        else
            write (ific, 120) lign1(1:80), lign1(81:160), lign1(161:nl11)
        end if
        if (nl22 .lt. 80) then
            write (ific, *) lign2(1:nl22)
        else if (nl22 .lt. 160) then
            write (ific, 116) lign2(1:80), lign2(81:nl22)
        else
            write (ific, 120) lign2(1:80), lign2(81:160), lign2(161:nl22)
        end if
!
        if (iret .eq. 0) then
        else if (iret .eq. 1) then
            call utmess('F', 'CALCULEL6_3')
        else if (iret .eq. 2) then
            call utmess('F', 'CALCULEL6_4')
        else if (iret .eq. 3) then
            call utmess('F', 'CALCULEL6_5')
        else
            call utmess('F', 'CALCULEL6_6')
        end if
        if (ctype(1:1) .ne. typr) then
            call utmess('F', 'CALCULEL6_8', sk=ctype)
        end if
!
        if (nk .ne. 0) then
!       cas des chaînes de caractères
            call tresu_str(tbtxt, zk80(irefk), valk, ific, .true._1)
            if (lref) then
                call tresu_str(tbtxt, zk80(irefkr), valk, ific, .false._1)
            end if
        else
!       cas des réels, entiers, complexes
            call tresu_print_all(tbtxt(1), tbtxt(2), .true._1, typr, nref, &
                                 crit, epsi, ssigne, zr(irefr), valr, &
                                 zi(irefi), vali, zc(irefc), valc, ignore=skip, &
                                 compare=ordgrd)
            if (lref) then
                call tresu_print_all(tbref(1), tbref(2), .false._1, typr, nref, &
                                     crit, epsir, ssigne, zr(irefrr), valr, &
                                     zi(irefir), vali, zc(irefcr), valc)
            end if
        end if
    end if
!
    if (nparfi .ne. 0) then
        call detrsd('TABLE', newta1)
    end if
    write (ific, *) ' '
!
100 format(/, 80('-'))
116 format(1x, a80, a)
120 format(1x, 2(a80), a)
!
    call jedema()
!
end subroutine
