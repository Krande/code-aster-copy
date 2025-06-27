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
subroutine trgene(ific, nocc)
!     COMMANDE:  TEST_RESU      MOT CLE FACTEUR "GENE"
! ----------------------------------------------------------------------
! aslint: disable=W1501

    use DynaGene_module

    implicit none
    integer(kind=8), intent(in) :: ific
    integer(kind=8), intent(in) :: nocc
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/codree.h"
#include "asterfort/extrac.h"
#include "asterfort/gettco.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tresu_ordgrd.h"
#include "asterfort/tresu_print_all.h"
#include "asterfort/tresu_read_refe.h"
#include "asterfort/trprec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zxtrac.h"
!
    character(len=6) :: nompro
    parameter(nompro='TRGENE')
    integer(kind=8) :: vali, iocc, iret, jlue, jordr, jdesc, jrefe, n1, n2, n3
    integer(kind=8) :: nbordr, numord, ncmp, nbinst, im, jcham, nbmode, jvecg
    integer(kind=8) :: jnume, jdeeq, istru, i, irefr, irefi, irefc, nref, nl1, nl2
    integer(kind=8) :: jfreq, nbfreq
    integer(kind=8) :: n1r, n2r, n3r, irefrr, irefir, irefcr, i_cham, i_bloc
    real(kind=8) :: valr, epsi, epsir, prec, temps, freq
    complex(kind=8) :: valc
    character(len=1) :: typres
    character(len=3) :: ssigne
    character(len=4) :: ch4
    character(len=8) :: crit, crit2, interp, mode
    character(len=11) :: motcle
    character(len=14) :: nugene
    character(len=16) :: nopara, nsym, k16b, tysd, ch16, tbtxt(2), tbref(2)
    character(len=19) :: cham19, knum, resu19
    character(len=24) :: travr, travi, travc, travrr, travir, travcr
    character(len=200) :: lign1, lign2
    aster_logical :: lref
    aster_logical :: skip
    real(kind=8) :: ordgrd
    integer(kind=8), pointer :: ordr(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    real(kind=8), pointer :: resu(:) => null()

    type(DynaGene) :: dyna_gene

!     ------------------------------------------------------------------
    call jemarq()
    motcle = 'RESU_GENE'
    travr = '&&'//nompro//'_TRAVR          '
    travi = '&&'//nompro//'_TRAVI          '
    travc = '&&'//nompro//'_TRAVC          '
    travrr = '&&'//nompro//'_TRAVR_R        '
    travir = '&&'//nompro//'_TRAVI_R        '
    travcr = '&&'//nompro//'_TRAVC_R        '
    irefi = 1
    irefr = 1
    irefc = 1
    irefir = 1
    irefcr = 1
    irefrr = 1
!
    do iocc = 1, nocc
        lign1 = ' '
        lign2 = ' '
!
        call trprec('GENE', iocc, epsi, crit, prec, &
                    crit2)
!
        call getvtx('GENE', 'VALE_ABS', iocc=iocc, scal=ssigne, nbret=n1)
!
        call getvr8('GENE', 'VALE_CALC', iocc=iocc, nbval=0, nbret=n1)
        call getvis('GENE', 'VALE_CALC_I', iocc=iocc, nbval=0, nbret=n2)
        call getvc8('GENE', 'VALE_CALC_C', iocc=iocc, nbval=0, nbret=n3)
        skip = .false.
        ordgrd = 1.d0
        if (n1 .ne. 0) then
            nref = -n1
            typres = 'R'
            call jedetr(travr)
            call wkvect(travr, 'V V R', nref, irefr)
            call getvr8('GENE', 'VALE_CALC', iocc=iocc, nbval=nref, vect=zr(irefr), &
                        nbret=iret)
            call tresu_ordgrd(zr(irefr), skip, ordgrd, mcf='GENE', iocc=iocc)
        else if (n2 .ne. 0) then
            nref = -n2
            typres = 'I'
            call jedetr(travi)
            call wkvect(travi, 'V V I', nref, irefi)
            call getvis('GENE', 'VALE_CALC_I', iocc=iocc, nbval=nref, vect=zi(irefi), &
                        nbret=iret)
        else if (n3 .ne. 0) then
            nref = -n3
            typres = 'C'
            call jedetr(travc)
            call wkvect(travc, 'V V C', nref, irefc)
            call getvc8('GENE', 'VALE_CALC_C', iocc=iocc, nbval=nref, vect=zc(irefc), &
                        nbret=iret)
        end if
! ----------------------------------------------------------------------
        lref = .false.
        call getvr8('GENE', 'PRECISION', iocc=iocc, scal=epsir, nbret=iret)
        if (iret .ne. 0) then
            lref = .true.
            call getvr8('GENE', 'VALE_REFE', iocc=iocc, nbval=0, nbret=n1r)
            call getvis('GENE', 'VALE_REFE_I', iocc=iocc, nbval=0, nbret=n2r)
            call getvc8('GENE', 'VALE_REFE_C', iocc=iocc, nbval=0, nbret=n3r)
            if (n1r .ne. 0) then
                ASSERT((n1r .eq. n1))
                nref = -n1r
                call jedetr(travrr)
                call wkvect(travrr, 'V V R', nref, irefrr)
                call getvr8('GENE', 'VALE_REFE', iocc=iocc, nbval=nref, vect=zr(irefrr), &
                            nbret=iret)
            else if (n2r .ne. 0) then
                ASSERT((n2r .eq. n2))
                nref = -n2r
                call jedetr(travir)
                call wkvect(travir, 'V V I', nref, irefir)
                call getvis('GENE', 'VALE_REFE_I', iocc=iocc, nbval=nref, vect=zi(irefir), &
                            nbret=iret)
            else if (n3r .ne. 0) then
                ASSERT((n3r .eq. n3))
                nref = -n3r
                call jedetr(travcr)
                call wkvect(travcr, 'V V C', nref, irefcr)
                call getvc8('GENE', 'VALE_REFE_C', iocc=iocc, nbval=nref, vect=zc(irefcr), &
                            nbret=iret)
            end if
        end if
        if (skip .and. .not. lref) then
            call utmess('I', 'TEST0_11')
        end if
! ----------------------------------------------------------------------
!
        call getvid('GENE', 'RESU_GENE', iocc=iocc, scal=resu19, nbret=n1)
        call gettco(resu19, tysd)
! ----------------------------------------------------------------------
        if (tysd .eq. 'VECT_ASSE_GENE') then
            call getvis('GENE', 'NUME_CMP_GENE', iocc=iocc, scal=ncmp, nbret=n1)
            call jeveuo(resu19//'.VALE', 'L', jlue)
            call jelira(resu19//'.VALE', 'TYPE', cval=k16b)
!
            call jeveuo(resu19//'.REFE', 'L', jrefe)
            mode = zk24(jrefe) (1:8)
            if (mode .eq. '        ') then
                nugene = zk24(jrefe+1) (1:14)
                call jeveuo(nugene//'.NUME.DEEQ', 'L', jdeeq)
                call jeveuo(nugene//'.NUME.NEQU', 'L', jnume)
                nbmode = zi(jnume)
                im = 0
                do i = 1, nbmode
                    istru = zi(jdeeq+2*(i-1)+2-1)
                    if (istru .lt. 0) goto 110
                    im = im+1
                    if (im .eq. ncmp) goto 114
110                 continue
                end do
                call utmess('F', 'CALCULEL6_98')
114             continue
                im = i
            else
                im = ncmp
            end if
!
            if (k16b(1:1) .ne. typres) then
                call utmess('F', 'CALCULEL6_95')
            else if (typres .eq. 'R') then
                valr = zr(jlue+im-1)
            else if (typres .eq. 'I') then
                vali = zi(jlue+im-1)
            else if (typres .eq. 'C') then
                valc = zc(jlue+im-1)
            end if
!
            lign1(1:21) = '---- '//motcle(1:9)
            lign1(22:22) = '.'
            lign2(1:21) = '     '//resu19(1:8)
            lign2(22:22) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_ORDRE'
            ch4 = '****'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_CMP_GENE'
            ch4 = ' '
            call codent(ncmp, 'G', ch4)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
!
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            if (nl1 .lt. 80) then
                write (ific, *) lign1(1:nl1)
            else if (nl1 .lt. 160) then
                write (ific, 1160) lign1(1:80), lign1(81:nl1)
            else
                write (ific, 1200) lign1(1:80), lign1(81:160), lign1( &
                    161:nl1)
            end if
            if (nl2 .lt. 80) then
                write (ific, *) lign2(1:nl2)
            else if (nl2 .lt. 160) then
                write (ific, 1160) lign2(1:80), lign2(81:nl2)
            else
                write (ific, 1200) lign2(1:80), lign2(81:160), lign2( &
                    161:nl2)
            end if
!
            call tresu_read_refe('GENE', iocc, tbtxt)
!
            if (lref) then
                tbref(1) = tbtxt(1)
                tbref(2) = tbtxt(2)
                tbtxt(1) = 'NON_REGRESSION'
            end if
!
            call tresu_print_all(tbtxt(1), tbtxt(2), .true._1, typres, nref, &
                                 crit, epsi, ssigne, zr(irefr), valr, &
                                 zi(irefi), vali, zc(irefc), valc, ignore=skip, &
                                 compare=ordgrd)
            if (lref) then
                call tresu_print_all(tbref(1), tbref(2), .false._1, typres, nref, &
                                     crit, epsir, ssigne, zr(irefrr), valr, &
                                     zi(irefir), vali, zc(irefcr), valc)
            end if
!
        else if (tysd .eq. 'MODE_GENE') then
!
            knum = '&&TRGENE.NUME_ORDRE'
            call rsutnu(resu19, 'GENE', iocc, knum, nbordr, &
                        prec, crit2, iret)
            if (iret .ne. 0) then
                call utmess('F', 'CALCULEL6_99', sk=resu19)
            end if
!
            call jeveuo(knum, 'L', jordr)
            numord = zi(jordr)
!
            call getvtx('GENE', 'PARA', iocc=iocc, scal=nopara, nbret=n1)
            if (n1 .ne. 0) then
                call rsadpa(resu19, 'L', 1, nopara, numord, &
                            1, sjv=jlue, styp=k16b)
                if (k16b(1:1) .ne. typres) then
                    call utmess('F', 'CALCULEL6_95')
                else if (typres .eq. 'R') then
                    valr = zr(jlue)
                else if (typres .eq. 'I') then
                    vali = zi(jlue)
                else if (typres .eq. 'C') then
                    valc = zc(jlue)
                end if
!
                lign1(1:21) = '---- '//motcle(1:9)
                lign1(22:22) = '.'
                lign2(1:21) = '     '//resu19(1:8)
                lign2(22:22) = '.'
                nl1 = lxlgut(lign1)
                nl2 = lxlgut(lign2)
                lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_ORDRE'
                ch4 = ' '
                call codent(numord, 'G', ch4)
                lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
                lign1(nl1+17:nl1+17) = '.'
                lign2(nl2+17:nl2+17) = '.'
                nl1 = lxlgut(lign1)
                nl2 = lxlgut(lign2)
                lign1(1:nl1+16) = lign1(1:nl1-1)//' PARA'
                lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nopara(1:16)
!
                nl1 = lxlgut(lign1)
                nl2 = lxlgut(lign2)
                if (nl1 .lt. 80) then
                    write (ific, *) lign1(1:nl1)
                else if (nl1 .lt. 160) then
                    write (ific, 1160) lign1(1:80), lign1(81:nl1)
                else
                    write (ific, 1200) lign1(1:80), lign1(81:160), &
                        lign1(161:nl1)
                end if
                if (nl2 .lt. 80) then
                    write (ific, *) lign2(1:nl2)
                else if (nl2 .lt. 160) then
                    write (ific, 1160) lign2(1:80), lign2(81:nl2)
                else
                    write (ific, 1200) lign2(1:80), lign2(81:160), &
                        lign2(161:nl2)
                end if
!
                call tresu_read_refe('GENE', iocc, tbtxt)
!
                if (lref) then
                    tbref(1) = tbtxt(1)
                    tbref(2) = tbtxt(2)
                    tbtxt(1) = 'NON_REGRESSION'
                end if
!
                call tresu_print_all(tbtxt(1), tbtxt(2), .true._1, typres, nref, &
                                     crit, epsi, ssigne, zr(irefr), valr, &
                                     zi(irefi), vali, zc(irefc), valc, ignore=skip, &
                                     compare=ordgrd)
                if (lref) then
                    call tresu_print_all(tbref(1), tbref(2), .false._1, typres, nref, &
                                         crit, epsir, ssigne, zr(irefrr), valr, &
                                         zi(irefir), vali, zc(irefcr), valc)
                end if
!
                call jedetr(knum)
                goto 100
            end if
!
            call getvtx('GENE', 'NOM_CHAM', iocc=iocc, scal=nsym, nbret=n1)
            call getvis('GENE', 'NUME_CMP_GENE', iocc=iocc, scal=ncmp, nbret=n1)
            call rsexch('F', resu19, nsym, numord, cham19, &
                        iret)
            call jeveuo(cham19//'.VALE_CALC', 'L', jlue)
            call jelira(cham19//'.VALE_CALC', 'TYPE', cval=k16b)
!
            call jeveuo(cham19//'.REFE', 'L', jrefe)
            mode = zk24(jrefe) (1:8)
            if (mode .eq. '        ') then
                nugene = zk24(jrefe+1) (1:14)
                call jeveuo(nugene//'.NUME.DEEQ', 'L', jdeeq)
                call jeveuo(nugene//'.NUME.NEQU', 'L', jnume)
                nbmode = zi(jnume)
                im = 0
                do i = 1, nbmode
                    istru = zi(jdeeq+2*(i-1)+2-1)
                    if (istru .lt. 0) goto 120
                    im = im+1
                    if (im .eq. ncmp) goto 124
120                 continue
                end do
                call utmess('F', 'CALCULEL6_98')
                goto 100
124             continue
                im = i
            else
                im = ncmp
            end if
!
            if (k16b(1:1) .ne. typres) then
                call utmess('F', 'CALCULEL6_95')
            else if (typres .eq. 'R') then
                valr = zr(jlue+im-1)
            else if (typres .eq. 'I') then
                vali = zi(jlue+im-1)
            else if (typres .eq. 'C') then
                valc = zc(jlue+im-1)
            end if
!
            lign1(1:21) = '---- '//motcle(1:9)
            lign1(22:22) = '.'
            lign2(1:21) = '     '//resu19(1:8)
            lign2(22:22) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_ORDRE'
            ch4 = ' '
            call codent(numord, 'G', ch4)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_CHAM'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nsym
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_CMP_GENE'
            ch4 = ' '
            call codent(ncmp, 'G', ch4)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
!
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            if (nl1 .lt. 80) then
                write (ific, *) lign1(1:nl1)
            else if (nl1 .lt. 160) then
                write (ific, 1160) lign1(1:80), lign1(81:nl1)
            else
                write (ific, 1200) lign1(1:80), lign1(81:160), lign1( &
                    161:nl1)
            end if
            if (nl2 .lt. 80) then
                write (ific, *) lign2(1:nl2)
            else if (nl2 .lt. 160) then
                write (ific, 1160) lign2(1:80), lign2(81:nl2)
            else
                write (ific, 1200) lign2(1:80), lign2(81:160), lign2( &
                    161:nl2)
            end if
!
            call tresu_read_refe('GENE', iocc, tbtxt)
!
            if (lref) then
                tbref(1) = tbtxt(1)
                tbref(2) = tbtxt(2)
                tbtxt(1) = 'NON_REGRESSION'
            end if
            call tresu_print_all(tbtxt(1), tbtxt(2), .true._1, typres, nref, &
                                 crit, epsi, ssigne, zr(irefr), valr, &
                                 zi(irefi), vali, zc(irefc), valc, ignore=skip, &
                                 compare=ordgrd)
            if (lref) then
                call tresu_print_all(tbref(1), tbref(2), .false._1, typres, nref, &
                                     crit, epsir, ssigne, zr(irefrr), valr, &
                                     zi(irefir), vali, zc(irefcr), valc)
            end if
            call jedetr(knum)
!
        else if (tysd .eq. 'HARM_GENE') then
            call getvtx('GENE', 'NOM_CHAM', iocc=iocc, scal=nsym, nbret=n1)
            call getvis('GENE', 'NUME_CMP_GENE', iocc=iocc, scal=ncmp, nbret=n1)
!
            interp = 'NON'
            call jeveuo(resu19//'.DISC', 'L', jfreq)
            call jelira(resu19//'.DISC', 'LONMAX', nbfreq)
!
            call getvr8('GENE', 'FREQ', iocc=iocc, scal=freq, nbret=n1)
            if (n1 .eq. 0) then
                call getvis('GENE', 'NUME_ORDRE', iocc=iocc, scal=numord, nbret=n1)
                call jeveuo(resu19//'.ORDR', 'L', jordr)
                do i = 1, nbfreq
                    if (numord .eq. zi(jordr+i-1)) then
                        freq = zr(jfreq+i-1)
                        exit
                    end if
                end do
            end if
!
            call jeexin(resu19//'.'//nsym(1:4), iret)
            if (iret .eq. 0) then
                call utmess('F', 'CALCULEL6_99', sk=resu19)
            end if
            call jeveuo(resu19//'.'//nsym(1:4), 'L', jcham)
            call jeveuo(resu19//'.DESC', 'L', jdesc)
            nbmode = zi(jdesc+2-1)
            call wkvect('&&TRGENE.CHAMP', 'V V C', nbmode, jvecg)
            call zxtrac(interp, prec, crit2, nbfreq, zr(jfreq), &
                        freq, zc(jcham), nbmode, zc(jvecg), iret)
            if (iret .ne. 0) then
                call utmess('F', 'CALCULEL6_2', sk=resu19)
            end if
            valc = zc(jvecg+ncmp-1)
!
            lign1(1:21) = '---- '//motcle(1:9)
            lign1(22:22) = '.'
            lign2(1:21) = '     '//resu19(1:8)
            lign2(22:22) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' FREQ'
            ch16 = ' '
            call codree(freq, 'E', ch16)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch16
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_CHAM'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nsym(1:4)
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_CMP_GENE'
            ch4 = ' '
            call codent(ncmp, 'G', ch4)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
!
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            if (nl1 .lt. 80) then
                write (ific, *) lign1(1:nl1)
            else if (nl1 .lt. 160) then
                write (ific, 1160) lign1(1:80), lign1(81:nl1)
            else
                write (ific, 1200) lign1(1:80), lign1(81:160), lign1( &
                    161:nl1)
            end if
            if (nl2 .lt. 80) then
                write (ific, *) lign2(1:nl2)
            else if (nl2 .lt. 160) then
                write (ific, 1160) lign2(1:80), lign2(81:nl2)
            else
                write (ific, 1200) lign2(1:80), lign2(81:160), lign2( &
                    161:nl2)
            end if
!
            call tresu_read_refe('GENE', iocc, tbtxt)
            if (lref) then
                tbref(1) = tbtxt(1)
                tbref(2) = tbtxt(2)
                tbtxt(1) = 'NON_REGRESSION'
            end if
            call tresu_print_all(tbtxt(1), tbtxt(2), .true._1, 'C', nref, &
                                 crit, epsi, ssigne, zr(irefr), valr, &
                                 zi(irefi), vali, zc(irefc), valc)
            if (lref) then
                call tresu_print_all(tbref(1), tbref(2), .false._1, 'C', nref, &
                                     crit, epsir, ssigne, zr(irefrr), valr, &
                                     zi(irefir), vali, zc(irefcr), valc)
            end if
            call jedetr('&&TRGENE.CHAMP')
!
        else if (tysd .eq. 'TRAN_GENE') then

            call dyna_gene%init(resu19(1:8))

            call getvtx('GENE', 'NOM_CHAM', iocc=iocc, scal=nsym, nbret=n1)
            call getvis('GENE', 'NUME_CMP_GENE', iocc=iocc, scal=ncmp, nbret=n1)
!
            interp = 'NON'
!
            call getvr8('GENE', 'INST', iocc=iocc, scal=temps, nbret=n1)
            if (n1 .eq. 0) then
                call getvis('GENE', 'NUME_ORDRE', iocc=iocc, scal=numord, nbret=n1)
                call dyna_gene%get_values_by_ordr(dyna_gene%ordr, numord, length=nbinst, vi=ordr)
                call dyna_gene%get_current_bloc(dyna_gene%ordr, i_bloc)
                call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)
                do i = 1, nbinst
                    if (numord .eq. ordr(i)) then
                        temps = disc(i)
                        exit
                    end if
                end do
            else
                call dyna_gene%get_values_by_disc(dyna_gene%disc, temps, length=nbinst, vr=disc)
                call dyna_gene%get_current_bloc(dyna_gene%disc, i_bloc)
            end if
!
            if (nsym(1:4) .eq. "DEPL") then
                i_cham = dyna_gene%depl
            else if (nsym(1:4) .eq. "VITE") then
                i_cham = dyna_gene%vite
            else
                ASSERT(.false.)
            end if
            call dyna_gene%get_current_bloc(dyna_gene%disc, i_bloc)
            call dyna_gene%get_values(i_cham, i_bloc, vr=resu)

            call jeveuo(resu19//'.DESC', 'L', jdesc)
            nbmode = zi(jdesc+2-1)
            call wkvect('&&TRGENE.CHAMP', 'V V R', nbmode, jvecg)
            call extrac(interp, prec, crit2, nbinst, disc, &
                        temps, resu, nbmode, zr(jvecg), iret)
            if (iret .ne. 0) then
                call utmess('F', 'CALCULEL6_2', sk=resu19)
            end if

            call dyna_gene%free

            valr = zr(jvecg+ncmp-1)
!
            lign1(1:21) = '---- '//motcle(1:9)
            lign1(22:22) = '.'
            lign2(1:21) = '     '//resu19(1:8)
            lign2(22:22) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' INST'
            ch16 = ' '
            call codree(temps, 'E', ch16)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch16
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_CHAM'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nsym(1:4)
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NUME_CMP_GENE'
            ch4 = ' '
            call codent(ncmp, 'G', ch4)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ch4
!
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            if (nl1 .lt. 80) then
                write (ific, *) lign1(1:nl1)
            else if (nl1 .lt. 160) then
                write (ific, 1160) lign1(1:80), lign1(81:nl1)
            else
                write (ific, 1200) lign1(1:80), lign1(81:160), lign1( &
                    161:nl1)
            end if
            if (nl2 .lt. 80) then
                write (ific, *) lign2(1:nl2)
            else if (nl2 .lt. 160) then
                write (ific, 1160) lign2(1:80), lign2(81:nl2)
            else
                write (ific, 1200) lign2(1:80), lign2(81:160), lign2( &
                    161:nl2)
            end if
!
            call tresu_read_refe('GENE', iocc, tbtxt)
            if (lref) then
                tbref(1) = tbtxt(1)
                tbref(2) = tbtxt(2)
                tbtxt(1) = 'NON_REGRESSION'
            end if
            call tresu_print_all(tbtxt(1), tbtxt(2), .true._1, 'R', nref, &
                                 crit, epsi, ssigne, zr(irefr), valr, &
                                 zi(irefi), vali, zc(irefc), valc, ignore=skip, &
                                 compare=ordgrd)
            if (lref) then
                call tresu_print_all(tbref(1), tbref(2), .false._1, 'R', nref, &
                                     crit, epsir, ssigne, zr(irefrr), valr, &
                                     zi(irefir), vali, zc(irefcr), valc)
            end if
            call jedetr('&&TRGENE.CHAMP')
        end if
        write (ific, *) ' '
100     continue
    end do
1160 format(1x, a80, a)
1200 format(1x, 2(a80), a)
    call jedetr(travr)
    call jedetr(travc)
    call jedetr(travi)
    call jedema()
end subroutine
