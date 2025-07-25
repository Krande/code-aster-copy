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

subroutine trchno(ific, nocc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tresu_champ_all.h"
#include "asterfort/tresu_champ_cmp.h"
#include "asterfort/tresu_champ_no.h"
#include "asterfort/tresu_ordgrd.h"
#include "asterfort/tresu_read_refe.h"
#include "asterfort/tresu_tole.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8), intent(in) :: ific
    integer(kind=8), intent(in) :: nocc
!     COMMANDE:  TEST_RESU
!                MOT CLE FACTEUR "CHAM_NO"
! ----------------------------------------------------------------------
!
    character(len=6) :: nompro
    parameter(nompro='TRCHNO')
!
    integer(kind=8) :: iocc, iret, nbcmp, n1, n2, n3, n4, ng
    integer(kind=8) :: n1r, n2r, n3r, irefrr, irefir, irefcr
    integer(kind=8) :: irefr, irefi, irefc, nref, nl1, nl2, nl11, nl22
    real(kind=8) :: epsi, epsir
    character(len=1) :: typres
    character(len=3) :: ssigne
    character(len=8) :: crit, noddl, nomma, typtes, exclgr
    character(len=11) :: motcle
    character(len=19) :: cham19
    character(len=16) :: tbtxt(2), tbref(2)
    character(len=33) :: nonoeu
    character(len=24) :: travr, travi, travc, travrr, travir, travcr, nogrno
    character(len=200) :: lign1, lign2
    aster_logical :: lref, l_parallel_mesh
    character(len=8), pointer :: nom_cmp(:) => null()
    aster_logical :: skip
    real(kind=8) :: ordgrd
!     ------------------------------------------------------------------
    call jemarq()
!
    motcle = 'CHAM_NO'
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
    irefrr = 1
    irefcr = 1
!
    do iocc = 1, nocc
        lign1 = ' '
        lign2 = ' '
        nonoeu = ' '
        noddl = ' '
        call getvid('CHAM_NO', 'CHAM_GD', iocc=iocc, scal=cham19, nbret=n1)
        lign1(1:21) = '---- '//motcle(1:8)
        lign1(22:22) = '.'
        lign2(1:21) = '     '//cham19(1:8)
        lign2(22:22) = '.'
!
        call tresu_read_refe('CHAM_NO', iocc, tbtxt)
!
        call getvtx('CHAM_NO', 'NOM_CMP', iocc=iocc, scal=noddl, nbret=n1)
        if (n1 .ne. 0) then
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_CMP'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//noddl
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
        end if
        call getvtx('CHAM_NO', 'VALE_ABS', iocc=iocc, scal=ssigne, nbret=n1)
        if (ssigne .eq. 'OUI') then
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' VALE_ABS'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//ssigne
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
        end if
        call tresu_tole(epsi, mcf='CHAM_NO', iocc=iocc)
        call getvtx('CHAM_NO', 'CRITERE', iocc=iocc, scal=crit, nbret=iret)
!
        call getvr8('CHAM_NO', 'VALE_CALC', iocc=iocc, nbval=0, nbret=n1)
        call getvis('CHAM_NO', 'VALE_CALC_I', iocc=iocc, nbval=0, nbret=n2)
        call getvc8('CHAM_NO', 'VALE_CALC_C', iocc=iocc, nbval=0, nbret=n3)
!
        skip = .false.
        ordgrd = 1.d0
        if (n1 .ne. 0) then
            nref = -n1
            typres = 'R'
            call jedetr(travr)
            call wkvect(travr, 'V V R', nref, irefr)
            call getvr8('CHAM_NO', 'VALE_CALC', iocc=iocc, nbval=nref, vect=zr(irefr), &
                        nbret=iret)
            call tresu_ordgrd(zr(irefr), skip, ordgrd, mcf='CHAM_NO', iocc=iocc)
        else if (n2 .ne. 0) then
            nref = -n2
            typres = 'I'
            call jedetr(travi)
            call wkvect(travi, 'V V I', nref, irefi)
            call getvis('CHAM_NO', 'VALE_CALC_I', iocc=iocc, nbval=nref, vect=zi(irefi), &
                        nbret=iret)
        else if (n3 .ne. 0) then
            nref = -n3
            typres = 'C'
            call jedetr(travc)
            call wkvect(travc, 'V V C', nref, irefc)
            call getvc8('CHAM_NO', 'VALE_CALC_C', iocc=iocc, nbval=nref, vect=zc(irefc), &
                        nbret=iret)
        end if
! ----------------------------------------------------------------------
        lref = .false.
        call getvr8('CHAM_NO', 'PRECISION', iocc=iocc, scal=epsir, nbret=iret)
        if (iret .ne. 0) then
            lref = .true.
            call getvr8('CHAM_NO', 'VALE_REFE', iocc=iocc, nbval=0, nbret=n1r)
            call getvis('CHAM_NO', 'VALE_REFE_I', iocc=iocc, nbval=0, nbret=n2r)
            call getvc8('CHAM_NO', 'VALE_REFE_C', iocc=iocc, nbval=0, nbret=n3r)
            if (n1r .ne. 0) then
                ASSERT((n1r .eq. n1))
                nref = -n1r
                call jedetr(travrr)
                call wkvect(travrr, 'V V R', nref, irefrr)
                call getvr8('CHAM_NO', 'VALE_REFE', iocc=iocc, nbval=nref, vect=zr(irefrr), &
                            nbret=iret)
            else if (n2r .ne. 0) then
                ASSERT((n2r .eq. n2))
                nref = -n2r
                call jedetr(travir)
                call wkvect(travir, 'V V I', nref, irefir)
                call getvis('CHAM_NO', 'VALE_REFE_I', iocc=iocc, nbval=nref, vect=zi(irefir), &
                            nbret=iret)
            else if (n3r .ne. 0) then
                ASSERT((n3r .eq. n3))
                nref = -n3r
                call jedetr(travcr)
                call wkvect(travcr, 'V V C', nref, irefcr)
                call getvc8('CHAM_NO', 'VALE_REFE_C', iocc=iocc, nbval=nref, vect=zc(irefcr), &
                            nbret=iret)
            end if
        end if
        if (skip .and. .not. lref) then
            call utmess('I', 'TEST0_11')
        end if
        if (skip .and. .not. lref) then
            call utmess('I', 'TEST0_11')
        end if
! ----------------------------------------------------------------------
!
        call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=nomma)
        l_parallel_mesh = isParallelMesh(nomma)
!
        call getvtx('CHAM_NO', 'TYPE_TEST', iocc=iocc, scal=typtes, nbret=n1)
!
        if (n1 .ne. 0) then
            !EXCLUS('NOEUD','GROUP_NO') avec 'TYPE_TEST'
            call getvtx('RESU', 'NOEUD', iocc=iocc, scal=exclgr, nbret=n2)
            if (n2 > 0 .and. l_parallel_mesh) then
                call utmess('F', 'MODELISA7_86')
            end if
            call getvtx('RESU', 'GROUP_NO', iocc=iocc, scal=exclgr, nbret=n3)
            if ((n2+n3) .gt. 0) then
                call utmess('A', 'CALCULEL6_96')
            end if

            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' TYPE_TEST'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//typtes
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
!
            call getvtx('CHAM_NO', 'NOM_CMP', iocc=iocc, nbval=0, nbret=n4)
            if (n4 .eq. 0) then
                nl1 = lxlgut(lign1)
                nl11 = lxlgut(lign1(1:nl1-1))
                nl2 = lxlgut(lign2)
                nl22 = lxlgut(lign2(1:nl2-1))
                if (nl11 .lt. 80) then
                    write (ific, *) lign1(1:nl11)
                else if (nl11 .lt. 160) then
                    write (ific, 116) lign1(1:80), lign1(81:nl11)
                else
                    write (ific, 120) lign1(1:80), lign1(81:160), &
                        lign1(161:nl11)
                end if
                if (nl22 .lt. 80) then
                    write (ific, *) lign2(1:nl22)
                else if (nl22 .lt. 160) then
                    write (ific, 116) lign2(1:80), lign2(81:nl22)
                else
                    write (ific, 120) lign2(1:80), lign2(81:160), &
                        lign2(161:nl22)
                end if
!
                if (lref) then
                    tbref(1) = tbtxt(1)
                    tbref(2) = tbtxt(2)
                    tbtxt(1) = 'NON_REGRESSION'
                end if
                call tresu_champ_all(cham19, typtes, typres, nref, tbtxt, &
                                     zi(irefi), zr(irefr), zc(irefc), epsi, crit, &
                                     .true._1, ssigne, ignore=skip, compare=ordgrd)
                if (lref) then
                    call tresu_champ_all(cham19, typtes, typres, nref, tbref, &
                                         zi(irefir), zr(irefrr), zc(irefcr), epsir, crit, &
                                         .false._1, ssigne)
                end if
!
            else
                nbcmp = -n4
                AS_ALLOCATE(vk8=nom_cmp, size=nbcmp)
                call getvtx('CHAM_NO', 'NOM_CMP', iocc=iocc, nbval=nbcmp, vect=nom_cmp, &
                            nbret=n4)
                if (lref) then
                    tbref(1) = tbtxt(1)
                    tbref(2) = tbtxt(2)
                    tbtxt(1) = 'NON_REGRESSION'
                end if
                call tresu_champ_cmp(cham19, typtes, typres, nref, tbtxt, &
                                     zi(irefi), zr(irefr), zc(irefc), epsi, lign1, &
                                     lign2, crit, ific, nbcmp, nom_cmp, &
                                     .true._1, ssigne, ignore=skip, compare=ordgrd)
                if (lref) then
                    call tresu_champ_cmp(cham19, typtes, typres, nref, tbref, &
                                         zi(irefir), zr(irefrr), zc(irefcr), epsir, lign1, &
                                         lign2, crit, ific, nbcmp, nom_cmp, &
                                         .false._1, ssigne)
                end if
                AS_DEALLOCATE(vk8=nom_cmp)
            end if
!
        else
!
            call getvtx('CHAM_NO', 'NOM_CMP', iocc=iocc, scal=noddl, nbret=n1)
            call getvem(nomma, 'NOEUD', 'CHAM_NO', 'NOEUD', iocc, &
                        1, nonoeu(1:8), n1)
            if (n1 .ne. 0) then
                nl1 = lxlgut(lign1)
                nl2 = lxlgut(lign2)
                lign1(1:nl1+16) = lign1(1:nl1-1)//' NOEUD'
                lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nonoeu(1:8)
!
                if (l_parallel_mesh) then
                    call utmess('F', 'MODELISA7_86')
                end if
            end if
!
            call getvem(nomma, 'GROUP_NO', 'CHAM_NO', 'GROUP_NO', iocc, &
                        1, nogrno, n2)

            ng = 0
            if (n2 == 0 .and. l_parallel_mesh) then
                call getvtx('RESU', 'GROUP_NO', iocc=iocc, nbval=1, scal=nogrno, nbret=ng)
            end if
            if ((n2 .ne. 0) .or. (ng .ne. 0)) then
                nl1 = lxlgut(lign1)
                nl2 = lxlgut(lign2)
                lign1(1:nl1+16) = lign1(1:nl1-1)//' GROUP_NO'
                lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nogrno
            end if
!
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            write (ific, *) lign1(1:nl1)
            write (ific, *) lign2(1:nl2)
!
!
            if (n1 .eq. 1) then
!            RIEN A FAIRE.
            else
                nonoeu = ' '
                if (n2 > 0) then
                    call utnono('F', nomma, 'NOEUD', nogrno, nonoeu(1:8), iret)
                end if
                nonoeu(10:33) = nogrno
            end if
!
            if (lref) then
                tbref(1) = tbtxt(1)
                tbref(2) = tbtxt(2)
                tbtxt(1) = 'NON_REGRESSION'
            end if
            call tresu_champ_no(cham19, nonoeu, noddl, nref, tbtxt, &
                                zi(irefi), zr(irefr), zc(irefc), typres, epsi, &
                                crit, .true._1, ssigne, ignore=skip, compare=ordgrd)
            if (lref) then
                call tresu_champ_no(cham19, nonoeu, noddl, nref, tbref, &
                                    zi(irefir), zr(irefrr), zc(irefcr), typres, epsir, &
                                    crit, .false._1, ssigne)
            end if
        end if
        write (ific, *) ' '
    end do
!
116 format(1x, a80, a)
120 format(1x, 2(a80), a)
!
    call jedetr(travr)
    call jedetr(travc)
    call jedetr(travi)
!
!
    call jedema()
end subroutine
