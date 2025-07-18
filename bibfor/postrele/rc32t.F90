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
subroutine rc32t()
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/compr8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rc32my.h"
#include "asterfort/rctres.h"
#include "asterfort/rcver1.h"
#include "asterfort/rcveri.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbliva.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
!
!     ------------------------------------------------------------------
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_ZE200 et B3200
!     AVEC LA METHODE DE SELECTION DES INSTANTS TRESCA
!     STOCKAGE DES CONTRAINTES TOTALES ET LINEARISEES :
!                 - THERMIQUES SOUS "RESU_THER"
!                 - DE PRESSION SOUS "RESU_PRES"
!                 - DUS AUX EFFORTS ET MOMENTS SOUS "RESU_MECA"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbsitu, nbther, nbpres, nbmeca, iocc
    character(len=24) :: jvorig, jvextr, valk(4)
    character(len=8) :: nocmp(6), crit(2), tabther, tabpres
    character(len=16) :: valek(2), methode
    integer(kind=8) :: nume1, n1, nume2, n2, nume3, n3, nn, ither, numether
    integer(kind=8) :: n5, ipres, numepres, imeca, numemeca, nbinst, jinst
    character(len=8) :: tabmeca, tableok, k8b, tabtemp
    aster_logical :: exist
    integer(kind=8) :: nbabsc, jabsc, ncmp, ndim, jorig, jextr, k, i, j
    real(kind=8) :: prec(2)
    real(kind=8) :: sigtot(1000*6), tminsn(2), tmaxsn(2), tminsp(2)
    real(kind=8) :: tmaxsp(2), tresc(2), trescb(2), r3(2), r3b(2)
    real(kind=8) :: tremin(2), treminb(2), tremax(2), tremaxb(2)
    real(kind=8) :: vale(2), momen0, momen1, siglin(1000*6), momen0ther, momen1ther
    real(kind=8) :: momen0pres, momen1pres, momen0mec, momen1mec, temp(8)
    real(kind=8) :: verifori, verifextr
    integer(kind=8) :: ibid, iret, kk, n4, jtemp, nb, l
    complex(kind=8) :: cbid
    real(kind=8), pointer :: contraintesth(:) => null()
    real(kind=8), pointer :: contraintespr(:) => null()
    real(kind=8), pointer :: contraintestot(:) => null()
    real(kind=8), pointer :: contraintesmec(:) => null()
    real(kind=8), pointer :: contraintesm(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    call getvtx(' ', 'METHODE', scal=methode, nbret=nb)
    if (methode .ne. 'TRESCA') goto 999
!
    call getfac('SITUATION', nbsitu)
!
    jvorig = '&&RC3200.TRANSIT.ORIG'
    jvextr = '&&RC3200.TRANSIT.EXTR'
    call jecrec(jvorig, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nbsitu)
    call jecrec(jvextr, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nbsitu)
    call jecrec('&&RC3200.TEMPCST', 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nbsitu)
!
    call getfac('RESU_THER', nbther)
    call getfac('RESU_PRES', nbpres)
    call getfac('RESU_MECA', nbmeca)
!
    nocmp(1) = 'SIXX'
    nocmp(2) = 'SIYY'
    nocmp(3) = 'SIZZ'
    nocmp(4) = 'SIXY'
    nocmp(5) = 'SIXZ'
    nocmp(6) = 'SIYZ'
!
    valek(1) = 'INST            '
    valek(2) = 'ABSC_CURV       '
!
    prec(1) = 1.0d-06
    prec(2) = 1.0d-06
    crit(1) = 'RELATIF'
    crit(2) = 'RELATIF'
!
    do i = 1, 8
        temp(i) = 0.d0
    end do
!
    l = 0
!
    do iocc = 1, nbsitu, 1
!
!------ on récupère les numéros des tables sous le mot clé situation
!------ puis les tables associées sous RESU_THER, RESU_PRES et RESU_MECA
        call getvis('SITUATION', 'NUME_RESU_THER', iocc=iocc, scal=nume1, nbret=n1)
        call getvis('SITUATION', 'NUME_RESU_PRES', iocc=iocc, scal=nume2, nbret=n2)
        call getvis('SITUATION', 'NUME_RESU_MECA', iocc=iocc, scal=nume3, nbret=n3)
        call getvid('SITUATION', 'TABL_TEMP', iocc=iocc, scal=tabtemp, nbret=n4)
        nn = n1+n2+n3
        if ((nn+n4) .eq. 0) goto 888
!
        if (n1 .ne. 0) then
            do ither = 1, nbther, 1
                call getvis('RESU_THER', 'NUME_RESU_THER', iocc=ither, scal=numether, nbret=n5)
                if (numether .eq. nume1) then
                    call getvid('RESU_THER', 'TABL_RESU_THER', iocc=ither, scal=tabther, &
                                nbret=n5)
                end if
            end do
        end if
!
        if (n2 .ne. 0) then
            do ipres = 1, nbpres, 1
                call getvis('RESU_PRES', 'NUME_RESU_PRES', iocc=ipres, scal=numepres, nbret=n5)
                if (numepres .eq. nume2) then
                    call getvid('RESU_PRES', 'TABL_RESU_PRES', iocc=ipres, scal=tabpres, &
                                nbret=n5)
                end if
            end do
        end if
!
        if (n3 .ne. 0) then
            do imeca = 1, nbmeca, 1
                call getvis('RESU_MECA', 'NUME_RESU_MECA', iocc=imeca, scal=numemeca, nbret=n5)
                if (numemeca .eq. nume3) then
                    call getvid('RESU_MECA', 'TABL_RESU_MECA', iocc=imeca, scal=tabmeca, &
                                nbret=n5)
                end if
            end do
        end if
!----------------------------------------------------
! ------ RECUPERATION DES INSTANTS ET DES ABSCISSES
!----------------------------------------------------
! ------ grace a la table thermique ou de pression ou mecanique
        if (n1 .ne. 0) then
            tableok = tabther
        else if (n2 .ne. 0) then
            tableok = tabpres
        else if (n3 .ne. 0) then
            tableok = tabmeca
        else
            tableok = tabtemp
        end if
!
! --------- on verifie l'ordre des noeuds de la table
        call rcveri(tableok)
!
! --------- on recupere les instants de la table
        call tbexip(tableok, valek(1), exist, k8b)
        if (.not. exist) then
            valk(1) = tableok
            valk(2) = valek(1)
            call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
        end if
        call tbexv1(tableok, valek(1), 'RC.INSTANT', 'V', nbinst, &
                    k8b)
        call jeveuo('RC.INSTANT', 'L', jinst)
!
! --------- on recupere les abscisses curvilignes de la table
        call tbexip(tableok, valek(2), exist, k8b)
        if (.not. exist) then
            valk(1) = tableok
            valk(2) = valek(2)
            call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
        end if
        call tbexv1(tableok, valek(2), 'RC.ABSC', 'V', nbabsc, &
                    k8b)
        call jeveuo('RC.ABSC', 'L', jabsc)
!
! --------- on vérifie que toutes les situations sont definies sur les mêmes abscisses
        l = l+1
        if (l .eq. 1) verifori = zr(jabsc)
        if (l .eq. 1) verifextr = zr(jabsc+nbabsc-1)
        if (l .gt. 1) then
            if (abs(verifori-zr(jabsc)) .gt. 1e-8) call utmess('F', 'POSTRCCM_45')
            if (abs(verifextr-zr(jabsc+nbabsc-1)) .gt. 1e-8) call utmess('F', 'POSTRCCM_45')
        end if
!
! ------ ON VERIFIE LA COHERENCE DES TABLES (THERMIQUE, PRESSION ET MECA)
        if (n1 .ne. 0) then
            if (n2 .ne. 0) then
                call rcver1('MECANIQUE', tabpres, tabther)
            end if
            if (n3 .ne. 0) then
                call rcver1('MECANIQUE', tabmeca, tabther)
            end if
        end if
        if (n2 .ne. 0) then
            if (n3 .ne. 0) then
                call rcver1('MECANIQUE', tabmeca, tabpres)
            end if
        end if
!
        if (nn .eq. 0) then
            if (n4 .ne. 0) then
                call jecroc(jexnum('&&RC3200.TEMPCST', iocc))
                call jeecra(jexnum('&&RC3200.TEMPCST', iocc), 'LONMAX', n4)
                call jeecra(jexnum('&&RC3200.TEMPCST', iocc), 'LONUTI', n4)
                call jeveuo(jexnum('&&RC3200.TEMPCST', iocc), 'E', jtemp)
!
                vale(1) = zr(jinst)
                vale(2) = zr(jabsc)
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, zr(jtemp), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            goto 888
        end if
!---------------------------------------------------------------
! ------ CREATION DES 2*8 TABLES A 6 composantes de contraintes
! ------ (totales et mécaniques en tminsp et tmaxsp)
! ------ (linéarisées totales + M_0 +flexion thermique en tminsn et tmaxsn)
! ------ (lin. ther. + lin. pres. + lin. meca en tminsn et tmaxsn)
! ------ tminsn et tmaxsn
!---------------------------------------------------------------
!
        AS_ALLOCATE(vr=contraintesth, size=nbabsc)
        AS_ALLOCATE(vr=contraintespr, size=nbabsc)
        AS_ALLOCATE(vr=contraintesmec, size=nbabsc)
        AS_ALLOCATE(vr=contraintestot, size=nbabsc)
        AS_ALLOCATE(vr=contraintesm, size=nbabsc)
!
        ncmp = 6
        ndim = 2*8*ncmp+4+4
        call jecroc(jexnum(jvorig, iocc))
        call jeecra(jexnum(jvorig, iocc), 'LONMAX', ndim)
        call jeecra(jexnum(jvorig, iocc), 'LONUTI', ndim)
        call jeveuo(jexnum(jvorig, iocc), 'E', jorig)
!
        call jecroc(jexnum(jvextr, iocc))
        call jeecra(jexnum(jvextr, iocc), 'LONMAX', ndim)
        call jeecra(jexnum(jvextr, iocc), 'LONUTI', ndim)
        call jeveuo(jexnum(jvextr, iocc), 'E', jextr)
!
        do k = 1, nbabsc*ncmp
            sigtot(k) = 0.d0
        end do
!
        do k = 1, 2
            tminsn(k) = 0.d0
            tmaxsn(k) = 0.d0
            tminsp(k) = 0.d0
            tmaxsp(k) = 0.d0
            tresc(k) = 0.d0
            trescb(k) = 0.d0
            r3(k) = 0.d0
            r3b(k) = 0.d0
            tremin(k) = 0.d0
            treminb(k) = 0.d0
            tremax(k) = 0.d0
            tremaxb(k) = 0.d0
        end do
!------------------------------------------------------------------
!        ETAPE 1 : Déterminer tminsp, tmaxsp, tminsn, tmasxn
!------------------------------------------------------------------
        do i = 1, nbinst
!
            vale(1) = zr(jinst+i-1)
!
            do j = 1, ncmp
!
                do k = 1, nbabsc
!
                    vale(2) = zr(jabsc+k-1)
!
                    if (n1 .ne. 0) then
                        call tbliva(tabther, 2, valek, [ibid], vale, &
                                    [cbid], k8b, crit, prec, nocmp(j), &
                                    k8b, ibid, contraintesth(k), cbid, k8b, &
                                    iret)
                        if (iret .ne. 0) then
                            valk(1) = tabther
                            valk(2) = nocmp(j)
                            valk(3) = valek(1)
                            valk(4) = valek(2)
                            call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                        valr=vale)
                        end if
                    end if
                    if (n2 .ne. 0) then
                        call tbliva(tabpres, 2, valek, [ibid], vale, &
                                    [cbid], k8b, crit, prec, nocmp(j), &
                                    k8b, ibid, contraintespr(k), cbid, k8b, &
                                    iret)
                        if (iret .ne. 0) then
                            valk(1) = tabpres
                            valk(2) = nocmp(j)
                            valk(3) = valek(1)
                            valk(4) = valek(2)
                            call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                        valr=vale)
                        end if
                    end if
                    if (n3 .ne. 0) then
                        call tbliva(tabmeca, 2, valek, [ibid], vale, &
                                    [cbid], k8b, crit, prec, nocmp(j), &
                                    k8b, ibid, contraintesmec(k), cbid, k8b, &
                                    iret)
                        if (iret .ne. 0) then
                            valk(1) = tabmeca
                            valk(2) = nocmp(j)
                            valk(3) = valek(1)
                            valk(4) = valek(2)
                            call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                        valr=vale)
                        end if
                    end if
                    contraintestot(k) = contraintesmec(k)+contraintespr(k)+contraintesth(k)
                    sigtot((k-1)*ncmp+j) = contraintesmec(k)+contraintespr(k)+contraintesth(k)
                end do
                call rc32my(nbabsc, zr(jabsc), contraintestot, momen0, momen1)
                siglin(1*j) = momen0-0.5d0*momen1
                siglin((nbabsc-1)*ncmp+j) = momen0+0.5d0*momen1
            end do
!
! --------- tmaxsp et tminsp
!
            kk = 0
            do k = 1, nbabsc, nbabsc-1
                kk = kk+1
                call rctres(sigtot((k-1)*ncmp+1), tresc(kk))
!
                if (trace(3, sigtot((k-1)*ncmp+1)) .lt. 0.d0) then
                    r3(kk) = tresc(kk)*(-1.0d0)
                else
                    r3(kk) = tresc(kk)
                end if
                if (i .eq. 1) then
                    tremin(kk) = r3(kk)
                    tremax(kk) = r3(kk)
                    tminsp(kk) = vale(1)
                    tmaxsp(kk) = vale(1)
                else
                    if (compr8(r3(kk), 'LT', tremin(kk), prec(1), 1)) then
                        tremin(kk) = r3(kk)
                        tminsp(kk) = vale(1)
                    end if
                    if (compr8(r3(kk), 'GT', tremax(kk), prec(1), 1)) then
                        tremax(kk) = r3(kk)
                        tmaxsp(kk) = vale(1)
                    end if
                end if
            end do
!
! --------- tmaxsn et tminsn
!
            kk = 0
!
            do k = 1, nbabsc, nbabsc-1
                kk = kk+1
                call rctres(siglin((k-1)*ncmp+1), trescb(kk))
!
                if (trace(3, siglin((k-1)*ncmp+1)) .lt. 0.d0) then
                    r3b(kk) = trescb(kk)*(-1.0d0)
                else
                    r3b(kk) = trescb(kk)
                end if
                if (i .eq. 1) then
                    treminb(kk) = r3b(kk)
                    tremaxb(kk) = r3b(kk)
                    tminsn(kk) = vale(1)
                    tmaxsn(kk) = vale(1)
                else
                    if (compr8(r3b(kk), 'LT', treminb(kk), prec(1), 1)) then
                        treminb(kk) = r3b(kk)
                        tminsn(kk) = vale(1)
                    end if
                    if (compr8(r3b(kk), 'GT', tremaxb(kk), prec(1), 1)) then
                        tremaxb(kk) = r3b(kk)
                        tmaxsn(kk) = vale(1)
                    end if
                end if
            end do
!
        end do
!
!------------------------------------------------------------------
!                 ETAPE 2 : Remplir les vecteurs
!------------------------------------------------------------------
! --------- pour tminsp a l'origine
        do j = 1, ncmp
            vale(1) = tminsp(1)
            vale(2) = zr(jabsc)
            if (n1 .ne. 0) then
                call tbliva(tabther, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesth(1), cbid, k8b, &
                            iret)
            end if
            if (n2 .ne. 0) then
                call tbliva(tabpres, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintespr(1), cbid, k8b, &
                            iret)
            end if
            if (n3 .ne. 0) then
                call tbliva(tabmeca, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesmec(1), cbid, k8b, &
                            iret)
            end if
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(1), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            contraintesm(1) = contraintespr(1)+contraintesmec(1)
            contraintestot(1) = contraintesth(1)+contraintespr(1)+contraintesmec(1)
!
            zr(jorig-1+j) = contraintestot(1)
            zr(jorig-1+72+j) = contraintesm(1)
            zr(jorig+90) = temp(1)
!
! --------- pour tmaxsp a l'origine
            vale(1) = tmaxsp(1)
            vale(2) = zr(jabsc)
            if (n1 .ne. 0) then
                call tbliva(tabther, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesth(1), cbid, k8b, &
                            iret)
            end if
            if (n2 .ne. 0) then
                call tbliva(tabpres, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintespr(1), cbid, k8b, &
                            iret)
            end if
            if (n3 .ne. 0) then
                call tbliva(tabmeca, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesmec(1), cbid, k8b, &
                            iret)
            end if
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(2), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            contraintesm(1) = contraintespr(1)+contraintesmec(1)
            contraintestot(1) = contraintesth(1)+contraintespr(1)+contraintesmec(1)
!
            zr(jorig-1+6+j) = contraintestot(1)
            zr(jorig-1+78+j) = contraintesm(1)
            zr(jorig+91) = temp(2)
!
! --------- pour tminsn a l'origine
            vale(1) = tminsn(1)
            do k = 1, nbabsc
                vale(2) = zr(jabsc+k-1)
                if (n1 .ne. 0) then
                    call tbliva(tabther, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesth(k), cbid, k8b, &
                                iret)
                end if
                if (n2 .ne. 0) then
                    call tbliva(tabpres, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintespr(k), cbid, k8b, &
                                iret)
                end if
                if (n3 .ne. 0) then
                    call tbliva(tabmeca, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesmec(k), cbid, k8b, &
                                iret)
                end if
                contraintestot(k) = contraintesth(k)+contraintespr(k)+contraintesmec(k)
            end do
!
            vale(2) = zr(jabsc)
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(3), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            zr(jorig+88) = temp(3)
!
            call rc32my(nbabsc, zr(jabsc), contraintestot, momen0, momen1)
            call rc32my(nbabsc, zr(jabsc), contraintesth, momen0ther, momen1ther)
            call rc32my(nbabsc, zr(jabsc), contraintespr, momen0pres, momen1pres)
            call rc32my(nbabsc, zr(jabsc), contraintesmec, momen0mec, momen1mec)
!
            zr(jorig-1+12+j) = momen0-0.5d0*momen1
            zr(jorig-1+24+j) = momen0pres
            zr(jorig-1+36+j) = momen0ther-0.5d0*momen1ther
            zr(jorig-1+48+j) = momen0pres-0.5d0*momen1pres
            zr(jorig-1+60+j) = momen0mec-0.5d0*momen1mec
            zr(jorig-1+92+j) = -0.5d0*momen1ther
!
! --------- pour tmaxsn a l'origine
            vale(1) = tmaxsn(1)
            do k = 1, nbabsc
                vale(2) = zr(jabsc+k-1)
!
                if (n1 .ne. 0) then
                    call tbliva(tabther, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesth(k), cbid, k8b, &
                                iret)
                end if
                if (n2 .ne. 0) then
                    call tbliva(tabpres, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintespr(k), cbid, k8b, &
                                iret)
                end if
                if (n3 .ne. 0) then
                    call tbliva(tabmeca, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesmec(k), cbid, k8b, &
                                iret)
                end if
                contraintestot(k) = contraintesth(k)+contraintespr(k)+contraintesmec(k)
            end do
!
            vale(2) = zr(jabsc)
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(4), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            zr(jorig+89) = temp(4)
!
            call rc32my(nbabsc, zr(jabsc), contraintestot, momen0, momen1)
            call rc32my(nbabsc, zr(jabsc), contraintesth, momen0ther, momen1ther)
            call rc32my(nbabsc, zr(jabsc), contraintespr, momen0pres, momen1pres)
            call rc32my(nbabsc, zr(jabsc), contraintesmec, momen0mec, momen1mec)
!
            zr(jorig-1+18+j) = momen0-0.5d0*momen1
            zr(jorig-1+30+j) = momen0pres
            zr(jorig-1+42+j) = momen0ther-0.5d0*momen1ther
            zr(jorig-1+54+j) = momen0pres-0.5d0*momen1pres
            zr(jorig-1+66+j) = momen0mec-0.5d0*momen1mec
            zr(jorig+84) = tminsn(1)
            zr(jorig+85) = tmaxsn(1)
            zr(jorig+86) = tminsp(1)
            zr(jorig+87) = tmaxsp(1)
            zr(jorig-1+98+j) = -0.5d0*momen1ther
!
! --------- pour tminsp a l'extremite
            vale(1) = tminsp(2)
            vale(2) = zr(jabsc+nbabsc-1)
            if (n1 .ne. 0) then
                call tbliva(tabther, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesth(nbabsc), cbid, k8b, &
                            iret)
            end if
            if (n2 .ne. 0) then
                call tbliva(tabpres, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintespr(nbabsc), cbid, k8b, &
                            iret)
            end if
            if (n3 .ne. 0) then
                call tbliva(tabmeca, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesmec(nbabsc), cbid, k8b, &
                            iret)
            end if
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(5), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            contraintesm(nbabsc) = contraintespr(nbabsc)+contraintesmec(nbabsc)
            contraintestot(nbabsc) = contraintesth(nbabsc)+contraintespr(nbabsc)+ &
                                     contraintesmec(nbabsc)
!
            zr(jextr-1+j) = contraintestot(nbabsc)
            zr(jextr-1+72+j) = contraintesm(nbabsc)
            zr(jextr+90) = temp(5)
!
! --------- pour tmaxsp a l'extremite
            vale(1) = tmaxsp(2)
            vale(2) = zr(jabsc+nbabsc-1)
            if (n1 .ne. 0) then
                call tbliva(tabther, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesth(nbabsc), cbid, k8b, &
                            iret)
            end if
            if (n2 .ne. 0) then
                call tbliva(tabpres, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintespr(nbabsc), cbid, k8b, &
                            iret)
            end if
            if (n3 .ne. 0) then
                call tbliva(tabmeca, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintesmec(nbabsc), cbid, k8b, &
                            iret)
            end if
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(6), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            contraintesm(nbabsc) = contraintespr(nbabsc)+contraintesmec(nbabsc)
            contraintestot(nbabsc) = contraintesth(nbabsc)+contraintespr(nbabsc)+ &
                                     contraintesmec(nbabsc)
!
            zr(jextr-1+6+j) = contraintestot(nbabsc)
            zr(jextr-1+78+j) = contraintesm(nbabsc)
            zr(jextr+91) = temp(6)
!
! --------- pour tminsn a l'extremite
            vale(1) = tminsn(2)
            do k = 1, nbabsc
                vale(2) = zr(jabsc+k-1)
                if (n1 .ne. 0) then
                    call tbliva(tabther, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesth(k), cbid, k8b, &
                                iret)
                end if
                if (n2 .ne. 0) then
                    call tbliva(tabpres, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintespr(k), cbid, k8b, &
                                iret)
                end if
                if (n3 .ne. 0) then
                    call tbliva(tabmeca, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesmec(k), cbid, k8b, &
                                iret)
                end if
                contraintestot(k) = contraintesth(k)+contraintespr(k)+contraintesmec(k)
            end do
            vale(2) = zr(jabsc+nbabsc-1)
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(7), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            zr(jextr+88) = temp(7)
!
            call rc32my(nbabsc, zr(jabsc), contraintestot, momen0, momen1)
            call rc32my(nbabsc, zr(jabsc), contraintesth, momen0ther, momen1ther)
            call rc32my(nbabsc, zr(jabsc), contraintespr, momen0pres, momen1pres)
            call rc32my(nbabsc, zr(jabsc), contraintesmec, momen0mec, momen1mec)
!
            zr(jextr-1+12+j) = momen0+0.5d0*momen1
            zr(jextr-1+24+j) = momen0pres
            zr(jextr-1+36+j) = momen0ther+0.5d0*momen1ther
            zr(jextr-1+48+j) = momen0pres+0.5d0*momen1pres
            zr(jextr-1+60+j) = momen0mec+0.5d0*momen1mec
            zr(jextr-1+92+j) = 0.5d0*momen1ther
!
! --------- pour tmaxsn a l'extremite
            vale(1) = tmaxsn(2)
            do k = 1, nbabsc
                vale(2) = zr(jabsc+k-1)
                if (n1 .ne. 0) then
                    call tbliva(tabther, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesth(k), cbid, k8b, &
                                iret)
                end if
                if (n2 .ne. 0) then
                    call tbliva(tabpres, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintespr(k), cbid, k8b, &
                                iret)
                end if
                if (n3 .ne. 0) then
                    call tbliva(tabmeca, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, contraintesmec(k), cbid, k8b, &
                                iret)
                end if
                contraintestot(k) = contraintesth(k)+contraintespr(k)+contraintesmec(k)
            end do
            vale(2) = zr(jabsc+nbabsc-1)
            if (n4 .ne. 0) then
                call tbliva(tabtemp, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, 'TEMP', &
                            k8b, ibid, temp(8), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tabtemp
                    valk(2) = 'TEMP'
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
            end if
            zr(jextr+89) = temp(8)
!
            call rc32my(nbabsc, zr(jabsc), contraintestot, momen0, momen1)
            call rc32my(nbabsc, zr(jabsc), contraintesth, momen0ther, momen1ther)
            call rc32my(nbabsc, zr(jabsc), contraintespr, momen0pres, momen1pres)
            call rc32my(nbabsc, zr(jabsc), contraintesmec, momen0mec, momen1mec)
!
            zr(jextr-1+18+j) = momen0+0.5d0*momen1
            zr(jextr-1+30+j) = momen0pres
            zr(jextr-1+42+j) = momen0ther+0.5d0*momen1ther
            zr(jextr-1+54+j) = momen0pres+0.5d0*momen1pres
            zr(jextr-1+66+j) = momen0mec+0.5d0*momen1mec
            zr(jextr+84) = tminsn(2)
            zr(jextr+85) = tmaxsn(2)
            zr(jextr+86) = tminsp(2)
            zr(jextr+87) = tmaxsp(2)
            zr(jextr-1+98+j) = 0.5d0*momen1ther
!
        end do
        call jedetr('RC.INSTANT')
        call jedetr('RC.ABSC')
        AS_DEALLOCATE(vr=contraintesm)
        AS_DEALLOCATE(vr=contraintesmec)
        AS_DEALLOCATE(vr=contraintesth)
        AS_DEALLOCATE(vr=contraintespr)
        AS_DEALLOCATE(vr=contraintestot)
!
888     continue
    end do
!
999 continue
    call jedema()
end subroutine
