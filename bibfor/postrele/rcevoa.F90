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
subroutine rcevoa(typtab, nommat)
!
! person_in_charge: the-hiep.chau at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rcvale.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nommat
    character(len=16) :: typtab
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     OPERATEUR POST_RCCM, TYPE_RESU_MECA = 'EVOLUTION'
!                                  OPTION = 'AMORCAGE'
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: n1, ibid, nbtran, nbpar, nbteta, nbcycl, nbabsc, nbins0, iret
    integer(kind=8) :: ind, i, k, l, jnock, jnocl, nk, nl, jteta, kinst, jfaij, jabsc
    integer(kind=8) :: i1, i2, ifm, niv, ndim, ioc, it, is1, is2, jinst, jtabl, jnbcy
    integer(kind=8) :: nbitot
    real(kind=8) :: r8b, prec(2), vale(4), valres(4), rij, rapp, fatot, fakl
    real(kind=8) :: fam, theta, rcal, sittef, aamorc, bamorc, damorc, ramorc, d
    real(kind=8) :: sitt1, sitt2, fkl
    complex(kind=8) :: cbid
    aster_logical :: exist, trouve
    integer(kind=8) :: icodre(4)
    character(len=8) :: k8b, nomres, crit(2), nompar, table
    character(len=16) :: motclf, valek(4), concep, nomcmd, nomval(4)
    character(len=19) :: nomf
    character(len=24) :: instan, ktheta, abscur, valk(7)
!
    integer(kind=8) :: nparm, npard
    parameter(nparm=2, npard=2)
    character(len=8) :: typarm(nparm), typard(npard)
    character(len=16) :: noparm(nparm), nopard(npard)
    data noparm/'THETA', 'FACT_AMORCAGE'/
    data typarm/'R', 'R'/
    data nopard/'THETA', 'FACT_AMORCAGE'/
    data typard/'R', 'R'/
! DEB ------------------------------------------------------------------
    call jemarq()
    r8b = 0.d0
!
    motclf = 'TRANSITOIRE'
    call getfac(motclf, nbtran)
    if (nbtran .eq. 0) goto 999
!
    call getres(nomres, concep, nomcmd)
!
    call infniv(ifm, niv)
!
    call tbcrsd(nomres, 'G')
    if (typtab .eq. 'VALE_MAX') then
        call tbajpa(nomres, nparm, noparm, typarm)
    else
        call tbajpa(nomres, npard, nopard, typard)
    end if
!
    valek(1) = 'ANGLE           '
    valek(2) = 'INST            '
    valek(3) = 'SIZZ            '
    valek(4) = 'ABSC_CURV       '
!
    prec(1) = 1.0d-06
    prec(2) = 1.0d-06
    crit(1) = 'RELATIF'
    crit(2) = 'RELATIF'
!
! --- RECUPERATION DES DONNEES MATERIAU
!
    nbpar = 0
    nompar = ' '
    nomval(1) = 'A_AMORC'
    nomval(2) = 'B_AMORC'
    nomval(3) = 'D_AMORC'
    nomval(4) = 'R_AMORC'
    call rcvale(nommat, 'RCCM', nbpar, nompar, [r8b], &
                4, nomval, valres, icodre, 2)
    aamorc = valres(1)
    bamorc = valres(2)
    damorc = valres(3)
    ramorc = valres(4)
!
    ktheta = '&&RCEVOA.THETA'
    instan = '&&RCEVOA.INSTANT'
    abscur = '&&RCEVOA.ABSC_CURV'
!
! --- LA PREMIERE TABLE DEFINIT LES THETA A TRAITER
!     ON VERIFIE QUE LES ABSC_CURV CORRESPONDENT AU RAMORC
!
    call getvid(motclf, 'TABL_SIGM_THETA', iocc=1, scal=table, nbret=n1)
    call tbexip(table, valek(1), exist, k8b)
    if (.not. exist) then
        valk(1) = table
        valk(2) = valek(1)
        call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
    end if
    call tbexv1(table, valek(1), ktheta, 'V', nbteta, &
                k8b)
    call jeveuo(ktheta, 'L', jteta)
!
    call tbexip(table, valek(4), exist, k8b)
    if (.not. exist) then
        valk(1) = table
        valk(2) = valek(4)
        call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
    end if
    call tbexv1(table, valek(4), abscur, 'V', nbabsc, &
                k8b)
    call jeveuo(abscur, 'L', jabsc)
!
! --- VERIFICATION DE LA DISTANCE D
!
    do it = 1, nbteta-1
        theta = (zr(jteta+it)-zr(jteta+it-1))*r8dgrd()
        d = zr(jabsc+it)-zr(jabsc+it-1)
        rcal = d/(2*sin(0.5d0*theta))
        rapp = abs((rcal-damorc)/damorc)
! ------ TOLERANCE DE 1%
        if (rapp .gt. 0.01d0) then
            vale(1) = rcal
            vale(2) = damorc
            call utmess('A', 'POSTRCCM_33', sk=table, nr=2, valr=vale)
        end if
    end do
!
! --- DETERMINATION DU NOMBRE DE SITUATION
!        = NOMBRE D'INSTANTS DES TRANSITOIRES
!
    call jecrec('&&RCEVOA.SITUATION', 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nbtran)
    nbitot = 0
    do ioc = 1, nbtran
!
        call getvid(motclf, 'TABL_SIGM_THETA', iocc=ioc, scal=table, nbret=n1)
        valk(1) = table
        do i1 = 1, 4
            call tbexip(table, valek(i1), exist, k8b)
            if (.not. exist) then
                valk(2) = valek(i1)
                call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
            end if
        end do
!
        call getvr8(motclf, 'INST', iocc=ioc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            nbins0 = -n1
            call jecroc(jexnum('&&RCEVOA.SITUATION', ioc))
            call jeecra(jexnum('&&RCEVOA.SITUATION', ioc), 'LONMAX', nbins0)
            call jeecra(jexnum('&&RCEVOA.SITUATION', ioc), 'LONUTI', nbins0)
            call jeveuo(jexnum('&&RCEVOA.SITUATION', ioc), 'E', kinst)
            call getvr8(motclf, 'INST', iocc=ioc, nbval=nbins0, vect=zr(kinst), &
                        nbret=n1)
        else
            call getvid(motclf, 'LIST_INST', iocc=ioc, scal=nomf, nbret=n1)
            if (n1 .ne. 0) then
                call jelira(nomf//'.VALE', 'LONMAX', nbins0)
                call jeveuo(nomf//'.VALE', 'L', jinst)
            else
                call tbexv1(table, valek(2), instan, 'V', nbins0, &
                            k8b)
                call jeveuo(instan, 'L', jinst)
            end if
            call jecroc(jexnum('&&RCEVOA.SITUATION', ioc))
            call jeecra(jexnum('&&RCEVOA.SITUATION', ioc), 'LONMAX', nbins0)
            call jeecra(jexnum('&&RCEVOA.SITUATION', ioc), 'LONUTI', nbins0)
            call jeveuo(jexnum('&&RCEVOA.SITUATION', ioc), 'E', kinst)
            do i = 1, nbins0
                zr(kinst-1+i) = zr(jinst-1+i)
            end do
            call jedetr(instan)
        end if
        nbitot = nbitot+nbins0
    end do
!
! --- CREATION DES OBJETS DE TRAVAIL
!
    call wkvect('&&RCEVOA.TABL_T', 'V V K8', nbitot, jtabl)
    call wkvect('&&RCEVOA.INST_T', 'V V R', nbitot, jinst)
    call wkvect('&&RCEVOA.NBCY_T', 'V V I', nbitot, jnbcy)
!
    ind = 0
    do ioc = 1, nbtran
!
        call jeveuo(jexnum('&&RCEVOA.SITUATION', ioc), 'L', kinst)
        call jelira(jexnum('&&RCEVOA.SITUATION', ioc), 'LONUTI', nbins0)
!
        call getvis(motclf, 'NB_OCCUR', iocc=ioc, scal=nbcycl, nbret=n1)
!
        call getvid(motclf, 'TABL_SIGM_THETA', iocc=ioc, scal=table, nbret=n1)
!
        do i = 1, nbins0
            ind = ind+1
            zk8(jtabl-1+ind) = table
            zr(jinst-1+ind) = zr(kinst-1+i)
            zi(jnbcy-1+ind) = nbcycl
        end do
    end do
!
! --- CALCUL DU FACTEUR D'AMORCAGE
!
    ndim = nbitot*nbitot
    call wkvect('&&RCEVFU.MATR_FA', 'V V R', ndim, jfaij)
    call wkvect('&&RCEVFU.NB_CYCL', 'V V I', nbitot, jnocl)
    call wkvect('&&RCEVFU.NB_CYCK', 'V V I', nbitot, jnock)
!
    do it = 1, nbteta
        vale(1) = zr(jteta+it-1)
        if (niv .eq. 2) then
            write (ifm, *) '   '
            write (ifm, *) '--->> ANGLE: ', zr(jteta+it-1)
        end if
!
        do i1 = 1, nbitot
!
            table = zk8(jtabl-1+i1)
            vale(2) = zr(jinst-1+i1)
            zi(jnock-1+i1) = zi(jnbcy-1+i1)
            zi(jnocl-1+i1) = zi(jnbcy-1+i1)
!
            call tbliva(table, 2, valek, [ibid], vale, &
                        [cbid], k8b, crit, prec, valek(3), &
                        k8b, ibid, sitt1, cbid, k8b, &
                        iret)
            if (iret .ne. 0) then
                valk(1) = table
                valk(2) = valek(3)
                valk(3) = valek(1)
                valk(4) = valek(2)
                call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                            valr=vale)
            end if
!
!  CORRECTION DE issue [23782] SUPPRESSION DES VALEURS ABOLUES sitt1 = abs(sitt1)
!
!
!
            do i2 = i1+1, nbitot
!
                table = zk8(jtabl-1+i2)
                vale(2) = zr(jinst-1+i2)
!
                call tbliva(table, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, valek(3), &
                            k8b, ibid, sitt2, cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = table
                    valk(2) = valek(3)
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
!
!  CORRECTION DE issue [23782]
!  SUPPRESSION DES VALEURS ABOLUES sitt2 = abs(sitt2)
!
                zr(jfaij-1+nbitot*(i1-1)+i2) = 0.d0
!
!CALCUL DU FACTEUR R
!EVOLUTION fiche [21804]  R = 0. si min(sitt1, sitt2) < 0.
!   on ne calcule R que si sitt1 et sitt2 sont positifs
!  d ou on n a plus besoin de valeur absolue ( issue [23782])
!
                if ((sitt1 .gt. r8prem()) .and. (sitt2 .gt. r8prem())) then
! ------------calcul du rapport de charge
                    rij = min(sitt1, sitt2)/max(sitt1, sitt2)
                else
                    rij = 0.
                end if
!
! ------------calcul de DELTASIGTT efficace
!
                sittef = abs(sitt1-sitt2)/(1.d0-(rij/ramorc))
!
! ------------calcul du facteur d amorcage elementaire
                fam = (sittef/aamorc)**(-1.d0/bamorc)
                zr(jfaij-1+nbitot*(i1-1)+i2) = fam
            end do
        end do
!
        fatot = 0.d0
!
        ind = 0
100     continue
        ind = ind+1
        if (niv .eq. 2) then
            if (ind .eq. 1) then
                write (ifm, *) 'MATRICE FACTEURS D''AMORCAGE INITIALE'
            else
                write (ifm, *) 'MATRICE FACTEURS D''AMORCAGE MODIFIEE'
            end if
            write (ifm, 1010) (zi(jnocl-1+l), l=1, nbitot)
            do k = 1, nbitot
                i1 = nbitot*(k-1)
                write (ifm, 1000) zi(jnock-1+k), (zr(jfaij-1+i1+l), l=1, &
                                                  nbitot)
            end do
        end if
!
        fam = 0.d0
        trouve = .false.
        do k = 1, nbitot
!
            if (zi(jnock-1+k) .eq. 0) goto 110
!
            do l = 1, nbitot
!
                if (zi(jnocl-1+l) .eq. 0) goto 112
!
                fakl = zr(jfaij-1+nbitot*(k-1)+l)
                if (fakl .gt. fam) then
                    trouve = .true.
                    fam = fakl
                    is1 = k
                    is2 = l
                    nl = zi(jnocl-1+l)
                    nk = zi(jnock-1+k)
                end if
!
112             continue
            end do
!
110         continue
        end do
!
        if (trouve) then
!
            nbcycl = min(nk, nl)
            fkl = fam*nbcycl
            if (niv .eq. 2) then
                write (ifm, 1020) '=> FACTEUR D''AMORCAGE MAXI: ', fam, &
                    is1, is2
                write (ifm, 1030) nbcycl, fkl
            end if
!
! -------- ON CUMULE
!
            fatot = fatot+fkl
!
! -------- ON MET A ZERO LES FACTEURS D'AMORCAGE INCRIMINES
!
            zi(jnocl-1+is2) = zi(jnocl-1+is2)-nbcycl
            zi(jnock-1+is1) = zi(jnock-1+is1)-nbcycl
            zi(jnocl-1+is1) = zi(jnocl-1+is1)-nbcycl
            zi(jnock-1+is2) = zi(jnock-1+is2)-nbcycl
            do i = 1, nbitot
                if (zi(jnock-1+is1) .eq. 0) then
                    zr(jfaij-1+nbitot*(is1-1)+i) = 0.d0
                end if
                if (zi(jnocl-1+is2) .eq. 0) then
                    zr(jfaij-1+nbitot*(i-1)+is2) = 0.d0
                end if
            end do
!
            goto 100
!
        end if
!
        if (niv .eq. 2) write (ifm, *) '-->> FACTEUR D''AMORCAGE CUMULE = ', fatot
!
        vale(2) = fatot
!
        if (typtab .eq. 'VALE_MAX') then
            call tbajli(nomres, nparm, noparm, [ibid], vale, &
                        [cbid], k8b, 0)
        else
            call tbajli(nomres, npard, nopard, [ibid], vale, &
                        [cbid], k8b, 0)
        end if
!
    end do
!
1000 format(1p, i10, '|', 40(e10.3, '|'))
1010 format(1p, ' NB_OCCUR ', '|', 40(i10, '|'))
1020 format(1p, a28, e12.5, ', LIGNE:', i4, ', COLONNE:', i4)
1030 format(1p, '   NB_OCCUR = ', i8, ', FA_KL = ', e9.2)
!
999 continue
!
    call jedema()
end subroutine
