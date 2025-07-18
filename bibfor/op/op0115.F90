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

subroutine op0115()
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8depi.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!     DEFINITION D'UNE MATRICE INTERSPECTRALE (DEFI_INTE_SPEC)
!     ------------------------------------------------------------------
    integer(kind=8) :: iocpf, iockt, ioccs
    integer(kind=8) :: ipf, ifonc, inum, ifreq, ikt, ics, ispec
    integer(kind=8) :: mxval, ibid, nbabs, nbfreq, nbval, nbfreqold
    integer(kind=8) :: lnumi, lnumj, lfonc, lfreq, lrefe, nbvalr
    integer(kind=8) :: lnoei, lnoej, lcmpi, lcmpj, n2, n3, n4, n5, n6, n7
!
    real(kind=8) :: valr, fmoy, ared, fmin, fmax, pas, freq, depi, num, den
    real(kind=8) :: rbid
!
    complex(kind=8) :: valc
    aster_logical :: diag
!
    character(len=8) :: nomu, fonc, k8bid, tfonc, nomref
    character(len=16) :: concep, nomcmd, motfac(3)
    character(len=19) :: chfonc
    character(len=24) :: chnumi, chnumj, chfreq, chvale
    character(len=24) :: chnoei, chnoej, chcmpi, chcmpj
    character(len=24), pointer :: prol(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call getres(nomu, concep, nomcmd)
!
    nomref = nomu(1:8)
!
    motfac(1) = 'PAR_FONCTION'
    motfac(2) = 'KANAI_TAJIMI'
    motfac(3) = 'CONSTANT'
!
    call getfac(motfac(1), iocpf)
    call getfac(motfac(2), iockt)
    call getfac(motfac(3), ioccs)
!
    mxval = iocpf+iockt+ioccs
!
    call wkvect(nomref//'.REFE', 'G V K16', 3, lrefe)
    zk16(lrefe) = 'DSP'
    zk16(lrefe+1) = 'TOUT'
    zk16(lrefe+2) = 'FREQ'
!
    call wkvect('&&FONC', 'V V K8', mxval, lfonc)
    chvale = nomref//'.VALE'
    call jecrec(chvale, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                mxval)
    chfreq = nomref//'.DISC'
!
    k8bid = 'BIDON'
    do ipf = 1, iocpf
        if (ipf .eq. 1) then
            call getvtx(motfac(1), 'NOEUD_I', iocc=ipf, nbval=0, nbret=n2)
            call getvis(motfac(1), 'NUME_ORDRE_I', iocc=ipf, nbval=0, nbret=n3)
!
            if (n2 .lt. 0) then
                chnoei = nomref//'.NOEI'
                chnoej = nomref//'.NOEJ'
                chcmpi = nomref//'.CMPI'
                chcmpj = nomref//'.CMPJ'
                call wkvect(chnoei, 'G V K8', mxval, lnoei)
                call wkvect(chnoej, 'G V K8', mxval, lnoej)
                call wkvect(chcmpi, 'G V K8', mxval, lcmpi)
                call wkvect(chcmpj, 'G V K8', mxval, lcmpj)
            else if (n3 .lt. 0) then
                chnumi = nomref//'.NUMI'
                chnumj = nomref//'.NUMJ'
                call wkvect(chnumi, 'G V I', mxval, lnumi)
                call wkvect(chnumj, 'G V I', mxval, lnumj)
            end if
        end if
!
        if (n2 .lt. 0) then
            call getvtx(motfac(1), 'NOEUD_I', iocc=ipf, scal=zk8(lnoei-1+ipf), nbret=nbval)
            call getvtx(motfac(1), 'NOEUD_J', iocc=ipf, scal=zk8(lnoej-1+ipf), nbret=nbval)
            call getvtx(motfac(1), 'NOM_CMP_I', iocc=ipf, scal=zk8(lcmpi-1+ipf), nbret=nbval)
            call getvtx(motfac(1), 'NOM_CMP_J', iocc=ipf, scal=zk8(lcmpj-1+ipf), nbret=nbval)
        else if (n3 .lt. 0) then
            call getvis(motfac(1), 'NUME_ORDRE_I', iocc=ipf, scal=zi(lnumi-1+ipf), nbret=nbval)
            call getvis(motfac(1), 'NUME_ORDRE_J', iocc=ipf, scal=zi(lnumj-1+ipf), nbret=nbval)
        end if
!
        call getvid(motfac(1), 'FONCTION', iocc=ipf, scal=zk8(lfonc-1+ipf), nbret=nbval)
    end do
!
    do ifonc = 1, iocpf
        fonc = zk8(lfonc-1+ifonc)
        chfonc = fonc//'           '
        call jeveuo(chfonc//'.VALE', 'L', vr=vale)
        call jelira(chfonc//'.VALE', 'LONMAX', nbval)
        call jeveuo(chfonc//'.PROL', 'L', vk24=prol)
        tfonc = prol(1) (1:8)
        if (tfonc .eq. 'FONCTION') nbfreq = nbval/2
        if (tfonc .eq. 'FONCT_C') nbfreq = nbval/3
        if (ifonc .eq. 1) then
            nbfreqold = nbfreq
        else
            if (nbfreq .ne. nbfreqold) then
                call utmess('F', 'ALGORITH3_79')
            end if
        end if
        diag = .false.
        if (n2 .lt. 0) then
            if ((zk8(lnoei-1+ifonc) .eq. zk8(lnoej-1+ifonc)) .and. &
                (zk8(lcmpi-1+ifonc) .eq. zk8(lcmpj-1+ifonc))) then
                nbabs = nbfreq
                diag = .true.
            else
                nbabs = nbfreq*2
            end if
        else if (n3 .lt. 0) then
            if (zi(lnumi-1+ifonc) .eq. zi(lnumj-1+ifonc)) then
                nbabs = nbfreq
                diag = .true.
            else
                nbabs = 2*nbfreq
            end if
        end if
        call jecroc(jexnum(chvale, ifonc))
        call jeecra(jexnum(chvale, ifonc), 'LONMAX', nbabs)
        call jeecra(jexnum(chvale, ifonc), 'LONUTI', nbabs)
        call jeveuo(jexnum(chvale, ifonc), 'E', ispec)
        if ((diag) .and. (tfonc .eq. 'FONCT_C')) then
            do inum = 1, nbabs
                zr(ispec-1+inum) = vale(nbfreq+2*(inum-1)+1)
            end do
        else if ((.not. diag) .and. (tfonc .eq. 'FONCTION')) then
            do inum = 1, nbfreq
                zr(ispec-1+2*inum-1) = vale(nbfreq+inum)
                zr(ispec-1+2*inum) = 0.d0
            end do
        else
            nbabs = nbval-nbfreq
            do inum = 1, nbabs
                zr(ispec-1+inum) = vale(nbfreq+inum)
            end do
        end if
    end do
!
    if (iocpf .gt. 0) then
        call jeexin(chfreq, ibid)
        if (ibid .eq. 0) then
            call wkvect(chfreq, 'G V R', nbfreq, lfreq)
            do ifreq = 1, nbfreq
                zr(lfreq-1+ifreq) = vale(ifreq)
            end do
        end if
    end if
!
    depi = r8depi()
!
    do ikt = 1, iockt
        if (ikt .eq. 1) then
            call getvtx(motfac(2), 'NOEUD_I', iocc=ikt, nbval=0, nbret=n4)
            call getvis(motfac(2), 'NUME_ORDRE_I', iocc=ikt, nbval=0, nbret=n5)
!
            if (n4 .lt. 0) then
                chnoei = nomref//'.NOEI'
                chnoej = nomref//'.NOEJ'
                chcmpi = nomref//'.CMPI'
                chcmpj = nomref//'.CMPJ'
                call jeexin(chnoei, n6)
                if (n6 .eq. 0) then
                    call wkvect(chnoei, 'G V K8', mxval, lnoei)
                    call wkvect(chnoej, 'G V K8', mxval, lnoej)
                    call wkvect(chcmpi, 'G V K8', mxval, lcmpi)
                    call wkvect(chcmpj, 'G V K8', mxval, lcmpj)
                end if
            else if (n5 .lt. 0) then
                chnumi = nomref//'.NUMI'
                chnumj = nomref//'.NUMJ'
                call jeexin(chnumi, n7)
                if (n7 .eq. 0) then
                    call wkvect(chnumi, 'G V I', mxval, lnumi)
                    call wkvect(chnumj, 'G V I', mxval, lnumj)
                end if
            end if
        end if
!
        if (n4 .lt. 0) then
            call getvtx(motfac(2), 'NOEUD_I', iocc=ikt, scal=zk8(lnoei-1+iocpf+ikt), &
                        nbret=nbval)
            call getvtx(motfac(2), 'NOEUD_J', iocc=ikt, scal=zk8(lnoej-1+iocpf+ikt), &
                        nbret=nbval)
            call getvtx(motfac(2), 'NOM_CMP_I', iocc=ikt, scal=zk8(lcmpi-1+iocpf+ikt), &
                        nbret=nbval)
            call getvtx(motfac(2), 'NOM_CMP_J', iocc=ikt, scal=zk8(lcmpj-1+iocpf+ikt), &
                        nbret=nbval)
        else if (n5 .lt. 0) then
            call getvis(motfac(2), 'NUME_ORDRE_I', iocc=ikt, scal=zi(lnumi-1+iocpf+ikt), &
                        nbret=nbval)
            call getvis(motfac(2), 'NUME_ORDRE_J', iocc=ikt, scal=zi(lnumj-1+iocpf+ikt), &
                        nbret=nbval)
        end if
        nbvalr = 0
        call getvr8(motfac(2), 'VALE_R', iocc=ikt, nbval=0, nbret=nbvalr)
        if (nbvalr .lt. 0) then
            call getvr8(motfac(2), 'VALE_R', iocc=ikt, scal=valr, nbret=nbval)
        else
            call getvc8(motfac(2), 'VALE_C', iocc=ikt, scal=valc, nbret=nbval)
! ON NE RETIENT QUE LA PARTIE REELLE
            valr = dble(valc)
        end if
        call getvr8(motfac(2), 'FREQ_MOY', iocc=ikt, scal=fmoy, nbret=nbval)
        call getvr8(motfac(2), 'AMOR_REDUIT', iocc=ikt, scal=ared, nbret=nbval)
        call getvr8(motfac(2), 'FREQ_MIN', iocc=ikt, scal=fmin, nbret=nbval)
        call getvr8(motfac(2), 'FREQ_MAX', iocc=ikt, scal=fmax, nbret=nbval)
        call getvr8(motfac(2), 'PAS', iocc=ikt, scal=pas, nbret=nbval)
        if (fmax .lt. fmin) then
            call utmess('F', 'SPECTRAL0_2', sk=motfac(2))
        end if
        nbfreq = int((fmax-fmin)/pas)+1
        ifonc = iocpf+ikt
        call jecroc(jexnum(chvale, ifonc))
        call jeecra(jexnum(chvale, ifonc), 'LONMAX', nbfreq)
        call jeecra(jexnum(chvale, ifonc), 'LONUTI', nbfreq)
        call jeveuo(jexnum(chvale, ifonc), 'E', ispec)
        do ifreq = 1, nbfreq
            freq = fmin+pas*(ifreq-1)
            if (ifreq .eq. nbfreq) freq = fmax
            rbid = 4.0d0*ared*ared*fmoy*fmoy*freq*freq
            num = rbid+fmoy*fmoy*fmoy*fmoy
            den = fmoy*fmoy-freq*freq
            den = den*den+rbid
            zr(ispec-1+ifreq) = depi*valr*num/den
        end do
    end do
!
    do ics = 1, ioccs
        if (ics .eq. 1) then
            call getvtx(motfac(3), 'NOEUD_I', iocc=ics, nbval=0, nbret=n6)
            call getvis(motfac(3), 'NUME_ORDRE_I', iocc=ics, nbval=0, nbret=n7)
!
            if (n6 .lt. 0) then
                chnoei = nomref//'.NOEI'
                chnoej = nomref//'.NOEJ'
                chcmpi = nomref//'.CMPI'
                chcmpj = nomref//'.CMPJ'
                call jeexin(chnoei, n4)
                if (n4 .eq. 0) then
                    call wkvect(chnoei, 'G V K8', mxval, lnoei)
                    call wkvect(chnoej, 'G V K8', mxval, lnoej)
                    call wkvect(chcmpi, 'G V K8', mxval, lcmpi)
                    call wkvect(chcmpj, 'G V K8', mxval, lcmpj)
                end if
            else if (n7 .lt. 0) then
                chnumi = nomref//'.NUMI'
                chnumj = nomref//'.NUMJ'
                call jeexin(chnumi, n5)
                if (n5 .eq. 0) then
                    call wkvect(chnumi, 'G V I', mxval, lnumi)
                    call wkvect(chnumj, 'G V I', mxval, lnumj)
                end if
            end if
        end if
!
        if (n6 .lt. 0) then
            call getvtx(motfac(3), 'NOEUD_I', iocc=ics, scal=zk8(lnoei-1+iocpf+iockt+ics), &
                        nbret=nbval)
            call getvtx(motfac(3), 'NOEUD_J', iocc=ics, scal=zk8(lnoej-1+iocpf+iockt+ics), &
                        nbret=nbval)
            call getvtx(motfac(3), 'NOM_CMP_I', iocc=ics, scal=zk8(lcmpi-1+iocpf+iockt+ics), &
                        nbret=nbval)
            call getvtx(motfac(3), 'NOM_CMP_J', iocc=ics, scal=zk8(lcmpj-1+iocpf+iockt+ics), &
                        nbret=nbval)
        else if (n7 .lt. 0) then
            call getvis(motfac(3), 'NUME_ORDRE_I', iocc=ics, scal=zi(lnumi-1+iocpf+iockt+ics), &
                        nbret=nbval)
            call getvis(motfac(3), 'NUME_ORDRE_J', iocc=ics, scal=zi(lnumj-1+iocpf+iockt+ics), &
                        nbret=nbval)
        end if
        ifonc = iocpf+iockt+ics
        nbvalr = 0
        call getvr8(motfac(3), 'VALE_R', iocc=ics, nbval=0, nbret=nbvalr)
        if (nbvalr .lt. 0) then
            call getvr8(motfac(3), 'VALE_R', iocc=ics, scal=valr, nbret=nbval)
        else
            call getvc8(motfac(3), 'VALE_C', iocc=ics, scal=valc, nbret=nbval)
        end if
        call getvr8(motfac(3), 'FREQ_MIN', iocc=ics, scal=fmin, nbret=nbval)
        call getvr8(motfac(3), 'FREQ_MAX', iocc=ics, scal=fmax, nbret=nbval)
        call getvr8(motfac(3), 'PAS', iocc=ics, scal=pas, nbret=nbval)
        if (fmax .lt. fmin) then
            call utmess('F', 'SPECTRAL0_2', sk=motfac(3))
        end if
        nbfreq = int((fmax-fmin)/pas)+1
        diag = .false.
        if (n6 .lt. 0) then
            if ((zk8(lnoei-1+ifonc) .eq. zk8(lnoej-1+ifonc)) .and. &
                (zk8(lcmpi-1+ifonc) .eq. zk8(lcmpj-1+ifonc))) then
                nbabs = nbfreq
                diag = .true.
            else
                nbabs = nbfreq*2
            end if
        else if (n7 .lt. 0) then
            if (zi(lnumi-1+ifonc) .eq. zi(lnumj-1+ifonc)) then
                nbabs = nbfreq
                diag = .true.
            else
                nbabs = nbfreq*2
            end if
        end if
        call jecroc(jexnum(chvale, ifonc))
        call jeecra(jexnum(chvale, ifonc), 'LONMAX', nbabs)
        call jeecra(jexnum(chvale, ifonc), 'LONUTI', nbabs)
        call jeveuo(jexnum(chvale, ifonc), 'E', ispec)
        do ifreq = 1, nbfreq
            if (diag) then
                if (nbvalr .lt. 0) then
                    zr(ispec-1+ifreq) = valr
                else
                    zr(ispec-1+ifreq) = dble(valc)
                end if
            else
                if (nbvalr .lt. 0) then
                    zr(ispec-1+2*ifreq-1) = valr
                    zr(ispec-1+2*ifreq) = 0.d0
                else
                    zr(ispec-1+2*ifreq-1) = dble(valc)
                    zr(ispec-1+2*ifreq) = dimag(valc)
                end if
            end if
        end do
    end do
!
    if ((iockt .gt. 0) .or. (ioccs .gt. 0)) then
        call jeexin(chfreq, ibid)
        if (ibid .eq. 0) then
            call wkvect(chfreq, 'G V R', nbfreq, lfreq)
            do ifreq = 1, nbfreq
                freq = fmin+pas*(ifreq-1)
                if (ifreq .eq. nbfreq) freq = fmax
                zr(lfreq-1+ifreq) = freq
            end do
        end if
    end if
!
    call titre()
!
    call jedema()
end subroutine
