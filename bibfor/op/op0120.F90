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
subroutine op0120()
    implicit none
!     CALCUL D'UNE MATRICE INTERSPECTRALE
!
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/calint.h"
#include "asterfort/fft.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/intimp.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rms.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifft, ifm, imatr, it, j, k
    integer(kind=8) :: kb, kf, kk, ktabl, l, l1, l2
    integer(kind=8) :: lcomp1, lfon, lint, long1, long2
    integer(kind=8) :: lresu1, lrms, ls1, lssx, ltabl, lvalc, lvale
    integer(kind=8) :: nbpts, nbpts2, nda, ndd
    integer(kind=8) :: nfcod, nfonc, niv, nmatr
    real(kind=8) :: bmatr, dfreq, dt, durana, durdec, frefin, freini
    real(kind=8) :: pts, pts1, pts2, pts3, tinst, tinst1
    real(kind=8) :: tinst2
!-----------------------------------------------------------------------
    integer(kind=8) :: long, ival(2), ier
    real(kind=8) :: resu, zero
    character(len=8) :: nomu, nomref
    character(len=16) :: concep, nomcmd
    character(len=19) :: nomfon
    character(len=24) :: nomobj
    character(len=24) :: chnumi, chnumj, chfreq, chvale
!
    integer(kind=8) :: ispec
    integer(kind=8) :: lnumi, lnumj, lfreq, lrefe, nbabs, mxval, ipf
!
!     ------------------------------------------------------------------
!
!     --- INITIALISATION DES DIVERS ---
    call jemarq()
!
    call getres(nomu, concep, nomcmd)
!
    call getvr8(' ', 'INST_INIT', scal=tinst1, nbret=l)
    call getvr8(' ', 'INST_FIN', scal=tinst2, nbret=l)
    call getvis(' ', 'NB_POIN', scal=nbpts, nbret=l)
    call getvid(' ', 'FONCTION', nbval=0, nbret=nfonc)
    nfonc = abs(nfonc)
!
!    --- VERIFICATION DU NOMBRE DE POINTS ---
    pts = log(dble(nbpts))/log(2.d0)
    pts1 = aint(pts)
    pts2 = abs(pts1-pts)
    pts3 = abs(1.d0-pts2)
    if (pts2 .ge. 1.d-06 .and. pts3 .ge. 1.d-06) then
        call utmess('F', 'ALGORITH9_56')
    end if
!
    call infmaj()
    call infniv(ifm, niv)
!
    nomref = nomu(1:8)
!
    call wkvect(nomref//'.REFE', 'G V K16', 3, lrefe)
    zk16(lrefe) = 'DSP'
    zk16(lrefe+1) = 'TOUT'
    zk16(lrefe+2) = 'FREQ'
!
    durana = tinst2-tinst1
    call getvr8(' ', 'DUREE_ANALYSE', scal=durana, nbret=nda)
!
    durdec = durana
    call getvr8(' ', 'DUREE_DECALAGE', scal=durdec, nbret=ndd)
!
    if (nda .ne. 0) then
        bmatr = ((tinst2-tinst1)-durana)/durdec
        nmatr = int(abs(bmatr)+1)
    else
        nmatr = 1
    end if
!
    call wkvect('&&OP0120.TEMP.LFON', 'V V K8', nfonc, lfon)
    call wkvect('&&OP0120.TEMP.VALE', 'V V C', nbpts, lvale)
!
    call getvid(' ', 'FONCTION', nbval=nfonc, vect=zk8(lfon), nbret=l)
!
    dt = durana/nbpts
    long = nbpts*nfonc/2
    nfcod = nfonc*(nfonc+1)/2
    long1 = nbpts*nfcod
    long2 = nmatr*nfcod
    nbpts2 = nbpts/2
    dfreq = 1.d0/durana
!C
    call wkvect('&&OP0120.TEMP.VALC', 'V V C', long, lvalc)
    call wkvect('&&OP0120.TEMP.LINT', 'V V R', nbpts, lint)
    call wkvect('&&OP0120.TEMP.LSSX', 'V V R', long1, lssx)
    call wkvect('&&OP0120.TEMP.LRMS', 'V V R', long2, lrms)
!C
    do imatr = 1, nmatr
        do kf = 1, nfonc
            nomfon = zk8(lfon+kf-1)
            do it = 1, nbpts
                tinst = tinst1+(imatr-1)*(durdec)+(it-1)*dt
                call fointe('F ', nomfon, 1, ['INST'], [tinst], &
                            resu, ier)
                zero = 0.d0
                zc(lvale+it-1) = dcmplx(resu, zero)
            end do
!       --- CALCUL DE LA TRANFORMEE DE FOURIER ---
            ifft = 1
            call fft(zc(lvale), nbpts, ifft)
!       ---
            do it = 1, nbpts2
                lresu1 = lvalc+(kf-1)*nbpts2+(it-1)
                zc(lresu1) = zc(lvale+it-1)*dt
            end do
        end do
        lcomp1 = 0
        do j = 1, nfonc
            do i = 1, j
!
                call calint(i, j, zc(lvalc), nbpts, zr(lint), &
                            long, durana)
!
                do kk = 1, nbpts2
                    l1 = lint+kk-1
                    l2 = lssx+kk-1+nbpts*lcomp1
                    zr(l2) = zr(l2)+zr(l1)
                    zr(l2+nbpts2) = zr(l2+nbpts2)+zr(l1+nbpts2)
                end do
                lcomp1 = lcomp1+1
            end do
        end do
        call rms(imatr, zr(lssx), long1, zr(lrms), long2, &
                 nbpts, nfcod, dfreq, nfonc)
    end do
    do kb = 1, long1
        ls1 = lssx+kb-1
        zr(ls1) = zr(ls1)/dble(nmatr)
    end do
!
!     --- CREATION DES NOMS DE FONCTIONS ---
    mxval = nfonc*(nfonc+1)/2
    chvale = nomref//'.VALE'
    call jecrec(chvale, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                mxval)
    chfreq = nomref//'.DISC'
    chnumi = nomref//'.NUMI'
    chnumj = nomref//'.NUMJ'
    call wkvect(chnumi, 'G V I', mxval, lnumi)
    call wkvect(chnumj, 'G V I', mxval, lnumj)
    call wkvect(chfreq, 'G V R', nbpts2, lfreq)
!
    do k = 1, nbpts2
        zr(lfreq+k-1) = (k-1)*dfreq
    end do
!
    ktabl = 1
    ipf = 0
    do j = 1, nfonc
        ival(2) = j
!
        do i = 1, j
            ival(1) = i
            ipf = ipf+1
            zi(lnumi-1+ipf) = ival(1)
            zi(lnumj-1+ipf) = ival(2)
!
            if (ival(1) .eq. ival(2)) then
                nbabs = nbpts2
            else
                nbabs = 2*nbpts2
            end if
!
            call jecroc(jexnum(chvale, ipf))
            call jeecra(jexnum(chvale, ipf), 'LONMAX', nbabs)
            call jeecra(jexnum(chvale, ipf), 'LONUTI', nbabs)
            call jeveuo(jexnum(chvale, ipf), 'E', ispec)
!
            do k = 1, nbpts2
                if (ival(1) .eq. ival(2)) then
                    l1 = ispec+k-1
                    l2 = lssx+nbpts*(ktabl-1)+k-1
                    zr(l1) = zr(l2)
                else
                    l1 = ispec+(k-1)*2
                    l2 = lssx+nbpts*(ktabl-1)+k-1
                    zr(l1) = zr(l2)
                    zr(l1+1) = zr(l2+nbpts2)
                end if
            end do
            ktabl = ktabl+1
        end do
    end do
!
    if (niv .ge. 1) then
        freini = 0.d0
        frefin = dfreq*(nbpts2-1)
        write (ifm, 200)
        write (ifm, 201) dfreq, freini, frefin
    end if
    if (niv .ge. 2) then
        nomobj = '&&OP0117.FONCTION'
        if (nfcod .ne. mxval) then
            call utmess('F', 'MODELISA2_89')
        end if
        call jeveuo(nomobj, 'L', ltabl)
        call intimp(ifm, zr(lrms), zk24(ltabl), nmatr, nfcod)
    end if
!
!
    call titre()
!
200 format('<PAS EN FREQUENCE>  <FREQ. INITIALE>  <FREQ. FINALE>')
201 format(4x, d11.4, 4x, d11.4, 4x, d11.4)
!
    call jedema()
end subroutine
