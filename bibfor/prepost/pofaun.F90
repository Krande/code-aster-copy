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
subroutine pofaun()
    implicit none
!     COMMANDE POST_FATIGUE :
!              CHARGEMENT PUREMENT UNIAXIAL
!     -----------------------------------------------------------------
!     ------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/fgcoke.h"
#include "asterfort/fgcorr.h"
#include "asterfort/fgcota.h"
#include "asterfort/fgdoba.h"
#include "asterfort/fgdohs.h"
#include "asterfort/fgdoma.h"
#include "asterfort/fgdomm.h"
#include "asterfort/fgdowh.h"
#include "asterfort/fgordo.h"
#include "asterfort/fgpeak.h"
#include "asterfort/fgpic2.h"
#include "asterfort/fgrain.h"
#include "asterfort/fgrccm.h"
#include "asterfort/fgrmax.h"
#include "asterfort/fgtahe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rccome.h"
#include "asterfort/rcpare.h"
#include "asterfort/rcvale.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbocc, ifonc, nbpts, i, n1, nbpapf, ifm, niv, nbp
    integer(kind=8) :: ivke, ivcorr, ivpoin, nbpoin, ivmax, ivmin, ivtrav
    integer(kind=8) :: ibid, intrav, ivpics, nbpics, nbcycl, nbpar, ivdome
    integer(kind=8) :: icodre(3), icodwo
    integer(kind=8) :: icodba, icodhs, icodma
    character(len=8) :: nomfon, result, txcum, k8b, nommat, kcorre
    character(len=8) :: method, nompar
    character(len=16) :: kdomm, nomcmd, methd1, nomres(3), cara
    character(len=32) :: pheno
    character(len=24) :: fvale
    real(kind=8) :: r8b, pseuil, rdomm, val(3), rampl
    complex(kind=8) :: cbid
    aster_logical :: lhaigh, fateps, lke
!     --- POST_FATI_UNIAX ----------------------------------------------
    parameter(nbpapf=5)
    character(len=1) :: typppf(nbpapf)
    character(len=16) :: nomppf(nbpapf)
    real(kind=8), pointer :: sigmax1(:) => null()
    real(kind=8), pointer :: sigmin1(:) => null()
    data nomppf/'CYCLE', 'VALE_MIN', 'VALE_MAX', 'DOMMAGE', 'DOMM_CUMU'/
    data typppf/'I', 'R', 'R', 'R', 'R'/
!     ------------------------------------------------------------------
!
    call jemarq()
!
    cbid = (0.d0, 0.d0)
    ibid = 0
    r8b = 0.d0
    fateps = .false.
    lhaigh = .false.
    lke = .false.
    ivpics = 0
    ivke = 0
    ivcorr = 0
!
!     --- RECUPERATION DU NIVEAU D'IMPRESSION ---
!
    call infniv(ifm, niv)
!
    call getres(result, k8b, nomcmd)
!
!     --- RECUPERATION DE LA FONCTION CHARGEMENT ---
!
    nomfon = ' '
    call getvid('HISTOIRE', 'SIGM', iocc=1, scal=nomfon, nbret=n1)
    call getvid('HISTOIRE', 'EPSI', iocc=1, scal=nomfon, nbret=n1)
    if (n1 .ne. 0) fateps = .true.
!
    fvale = nomfon//'           .VALE'
    call jelira(fvale, 'LONMAX', nbpts)
    call jeveuo(fvale, 'L', ifonc)
    nbpts = nbpts/2
    call wkvect('&&POFAUN.FONC.POIN', 'V V R', nbpts, ivpoin)
!
!     --- IMPRESSION DE LA FONCTION ----
    if (niv .eq. 2) then
        write (ifm, '(1X,A)') 'VALEURS DE LA FONCTION CHARGEMENT:'
        do i = 1, nbpts
            write (ifm, 1000) zr(ifonc+i-1), zr(ifonc+nbpts+i-1)
        end do
    end if
!
!
!     --- RECUPERATION DU COEFFICIENT D'AMPLIFICATION ---
    rampl = 1
    call getfac('COEF_MULT', nbocc)
    if (nbocc .ne. 0) then
        call getvr8('COEF_MULT', 'KT', iocc=1, scal=rampl, nbret=n1)
!        CALL FGAMPL(RAMPL,NBPTS,ZR(NBPTS+IFONC))
!
    end if
!
!     --- EXTRACTION DES PICS DE LA FONCTION DE CHARGEMENT ---
!
    call getvr8(' ', 'DELTA_OSCI', scal=pseuil, nbret=n1)
    call fgpeak(nomfon, pseuil, rampl, nbpoin, zr(ivpoin))
!
!     --- IMPRESSION DES PICS EXTRAITS DE LA FONCTION ----
    if (niv .eq. 2) then
        write (ifm, *)
        write (ifm, '(1X,A)') 'PICS EXTRAITS DE LA FONCTION CHARGEMENT'
        write (ifm, '(1X,A)') 'APRES AVOIR PRIS EN COMPTE DE KT'
        write (ifm, 1010) pseuil, nbpoin
        write (ifm, *)
        write (ifm, '(4(1X,E18.6))') (zr(ivpoin+i-1), i=1, nbpoin)
    end if
!
!
!     ---RECUPERATION DE LA LOI DE COMPTAGES DE CYCLES
!
    call getvtx(' ', 'COMPTAGE', scal=methd1, nbret=n1)
    if (methd1(9:12) .ne. '_MAX') then
        method = methd1(1:8)
    else
        method = 'RFLO_MAX'
    end if
!
    call wkvect('&&POFAUN.SIGMAX', 'V V R', nbpoin+2, ivmax)
    call wkvect('&&POFAUN.SIGMIN', 'V V R', nbpoin+2, ivmin)
    AS_ALLOCATE(vr=sigmax1, size=nbpoin+2)
    AS_ALLOCATE(vr=sigmin1, size=nbpoin+2)
    call wkvect('&&POFAUN.POIN.TRAV', 'V V R', nbpoin+2, ivtrav)
    call wkvect('&&POFAUN.NUME.TRAV', 'V V I', 2*(nbpoin+2), intrav)
    if (method .eq. 'RAINFLOW') then
        call wkvect('&&POFAUN.FONC.PICS', 'V V R', nbpoin+2, ivpics)
        call fgpic2(method, zr(ivtrav), zr(ivpoin), nbpoin, zr(ivpics), &
                    nbpics)
        call fgrain(zr(ivpics), nbpics, zi(intrav), nbcycl, zr(ivmin), &
                    zr(ivmax))
    else if (method .eq. 'RFLO_MAX') then
!
        call wkvect('&&POFAUN.FONC.PICS', 'V V R', nbpoin+2, ivpics)
        call fgpic2(method, zr(ivtrav), zr(ivpoin), nbpoin, zr(ivpics), &
                    nbpics)
        call fgrain(zr(ivpics), nbpics, zi(intrav), nbcycl, sigmin1, &
                    sigmax1)
!
        call fgrmax(nbcycl, sigmin1, sigmax1, zr(ivmin), zr(ivmax))
!
    else if (method .eq. 'RCCM') then
        call fgordo(nbpoin, zr(ivpoin), zr(ivtrav))
        call fgrccm(nbpoin, zr(ivtrav), nbcycl, zr(ivmin), zr(ivmax))
    else if (method .eq. 'NATUREL') then
        call fgcota(nbpoin, zr(ivpoin), nbcycl, zr(ivmin), zr(ivmax))
    else
        call utmess('F', 'FATIGUE1_15')
    end if
    if (nbcycl .eq. 0) then
        call utmess('F', 'FATIGUE1_16')
    end if
!
!     --- CORRECTION ELASTO-PLASTIQUE ---
!
    kcorre = ' '
    call getvtx(' ', 'CORR_KE', scal=kcorre, nbret=n1)
    call getvid(' ', 'MATER', scal=nommat, nbret=n1)
    if (kcorre .eq. 'RCCM') then
        nomres(1) = 'N_KE'
        nomres(2) = 'M_KE'
        nomres(3) = 'SM'
        nbpar = 0
        nompar = ' '
        call rcvale(nommat, 'RCCM', nbpar, nompar, [r8b], &
                    3, nomres, val, icodre, 2)
        call wkvect('&&POFAUN.KE', 'V V R', nbcycl, ivke)
        lke = .true.
        call fgcoke(nbcycl, zr(ivmin), zr(ivmax), val(1), val(2), &
                    val(3), zr(ivke))
    end if
!
!     --- CALCUL DU DOMMAGE ELEMENTAIRE ---
!
    kdomm = ' '
    call getvtx(' ', 'DOMMAGE', scal=kdomm, nbret=n1)
!
    call wkvect('&&POFAUN.DOMM.ELEM', 'V V R', nbcycl, ivdome)
!
!     --- CALCUL DU DOMMAGE ELEMENTAIRE DE WOHLER ---
!         ---------------------------------------
    if (kdomm .eq. 'WOHLER') then
!        ---CORRECTION DE HAIG (GOODMANN OU GERBER)
        kcorre = ' '
        call getvtx(' ', 'CORR_SIGM_MOYE', scal=kcorre, nbret=n1)
        if (kcorre .ne. ' ') then
            nomres(1) = 'SU'
            nbpar = 0
            nompar = ' '
            call rcvale(nommat, 'RCCM', nbpar, nompar, [r8b], &
                        1, nomres, val, icodre, 2)
            call wkvect('&&POFAUN.HAIG', 'V V R', nbcycl, ivcorr)
            lhaigh = .true.
            call fgcorr(nbcycl, zr(ivmin), zr(ivmax), kcorre, val(1), &
                        zr(ivcorr))
        end if
!
        pheno = 'FATIGUE'
        call rccome(nommat, pheno, icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_24')
        end if
        cara = 'WOHLER'
        call rcpare(nommat, pheno, cara, icodwo)
        cara = 'A_BASQUIN'
        call rcpare(nommat, pheno, cara, icodba)
        cara = 'A0'
        call rcpare(nommat, pheno, cara, icodhs)
        if (icodwo .eq. 0) then
            call fgdowh(nommat, nbcycl, zr(ivmin), zr(ivmax), lke, &
                        zr(ivke), lhaigh, zr(ivcorr), zr(ivdome))
        else if (icodba .eq. 0) then
            call fgdoba(nommat, nbcycl, zr(ivmin), zr(ivmax), lke, &
                        zr(ivke), lhaigh, zr(ivcorr), zr(ivdome))
        else if (icodhs .eq. 0) then
            call fgdohs(nommat, nbcycl, zr(ivmin), zr(ivmax), lke, &
                        zr(ivke), lhaigh, zr(ivcorr), zr(ivdome))
        end if
!
!     --- CALCUL DU DOMMAGE ELEMENTAIRE DE MANSON_COFFIN ----
!         ----------------------------------------------
    else if (kdomm .eq. 'MANSON_COFFIN') then
        if (.not. fateps) then
            call utmess('F', 'FATIGUE1_17')
        end if
        pheno = 'FATIGUE'
        call rccome(nommat, pheno, icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_24')
        end if
        cara = 'MANSON_COFFIN'
        call rcpare(nommat, pheno, cara, icodma)
        if (icodma .eq. 0) then
            call fgdoma(nommat, nbcycl, zr(ivmin), zr(ivmax), zr(ivdome))
        else
            call utmess('F', 'FATIGUE1_18')
        end if
!
!     --- CALCUL DU DOMMAGE ELEMENTAIRE DE TAHERI ---
!         ---------------------------------------
    else if (kdomm(1:6) .eq. 'TAHERI') then
        if (fateps) then
            call fgtahe(kdomm, nbcycl, zr(ivmin), zr(ivmax), zr(ivdome))
        else
            call utmess('F', 'FATIGUE1_19')
        end if
!
    else if (kdomm .eq. ' ') then
    else
        call utmess('F', 'FATIGUE1_20')
    end if
!
!     --- CREATION DE LA TABLE ---
!
    call tbcrsd(result, 'G')
    call tbajpa(result, nbpapf, nomppf, typppf)
!
    nbp = 4
    if (kdomm .eq. ' ') nbp = 3
    do i = 1, nbcycl
        val(1) = zr(ivmin+i-1)
        val(2) = zr(ivmax+i-1)
        val(3) = zr(ivdome+i-1)
        call tbajli(result, nbp, nomppf, [i], val, &
                    [cbid], k8b, 0)
    end do
!
!     --- CALCUL DU DOMMAGE TOTAL ---
!
    txcum = ' '
    call getvtx(' ', 'CUMUL', scal=txcum, nbret=n1)
    if (txcum .eq. 'LINEAIRE') then
!
        call fgdomm(nbcycl, zr(ivdome), rdomm)
!
        call tbajli(result, 1, nomppf(5), [ibid], [rdomm], &
                    [cbid], k8b, 0)
!
    end if
!
    call jedetr('&&POFAUN.FONC.POIN')
    call jedetr('&&POFAUN.SIGMAX')
    call jedetr('&&POFAUN.SIGMIN')
    AS_DEALLOCATE(vr=sigmax1)
    AS_DEALLOCATE(vr=sigmin1)
    call jedetr('&&POFAUN.POIN.TRAV')
    call jedetr('&&POFAUN.NUME.TRAV')
    call jedetr('&&POFAUN.DOMM.ELEM')
    if (ivpics .ne. 0) call jedetr('&&POFAUN.FONC.PICS')
    if (ivke .ne. 0) call jedetr('&&POFAUN.KE')
    if (ivcorr .ne. 0) call jedetr('&&POFAUN.HAIG')
!
1000 format(2x, e18.6, 5x, e18.6)
1010 format(1x, 'SEUIL = ', e18.6, 10x, 'NB DE PICS = ', i5)
!
    call jedema()
end subroutine
