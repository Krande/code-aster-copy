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
subroutine rvechn(ssch19, sdlieu, sdeval)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/tremno.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=19) :: ssch19, sdlieu, sdeval
!
!     OPERATION D' EXTRACTION DU POST-TRAITEMENT SUR UNE LISTE DE NOEUDS
!     ------------------------------------------------------------------
! IN  SSCH19 : K : NOM DU SOUS CHAMP DE GRANDEUR
! IN  SDLIEU : K : NOM DE LA SD REPRESENTANT LE LIEU
! OUT SDEVAL : K : NOM DE LA SD SOUS_CHAMP_GD PRODUITES
!            :   :(DESCRIPTION : CF RVPSTE)
!     ------------------------------------------------------------------
!
!
!
    character(len=24) :: invale, inpadr, inpcmp, innoma, innugd, inpnco, inpnsp
    character(len=24) :: ouvale, oupadr, oupcmp, ouerre, ounoma, oupnbn, ounugd
    character(len=24) :: oupnb2
    character(len=24) :: nrefe, nabsc, ndesc, nnumnd, nindir, nnmail, oupnco
    character(len=24) :: oupnsp
    character(len=24) :: valk
    character(len=19) :: sdemno
    character(len=15) :: nrepma
    character(len=8) :: mailla
    character(len=4) :: docu
!
    integer(kind=8) :: aipadr, aipcmp, iocer, j, anmail, anuma, aipnco, aipnsp, aopnco
    integer(kind=8) :: aopnbn, aovale, aopadr, aopcmp, aoerre, ainugd, aivale, aopnsp
    integer(kind=8) :: arefe, adesc, nbcmp, i, anumnd, acmpgd, lpt, aopnb2
    integer(kind=8) :: nbmpst, nbnpst, nbocer, n, m, adrin, adrou, nbm, numm
    integer(kind=8) :: nbtcmp, sdvacp, aindir, pt, nsp, nco, lmc, lcc, lsc, lms
    integer(kind=8) :: indi1, indi2, ier
    integer(kind=8) :: vali, ilong, k, l, lnc, ncom, nspm
!
    aster_logical :: trouve, lnomnoe, lnommai
!
    character(len=1) :: cbid
    integer(kind=8), pointer :: nund(:) => null()
    data cbid/' '/
!
!==================== CORPS DE LA ROUTINE =============================
!
    call jemarq()
    nnumnd = '&&RVECHN.NUM.NOEUD.LISTE'
    sdemno = '&&RVECHN.SDEMNO   '
    invale = ssch19//'.VALE'
    inpadr = ssch19//'.PADR'
    inpcmp = ssch19//'.PCMP'
    inpnco = ssch19//'.PNCO'
    inpnsp = ssch19//'.PNSP'
    innoma = ssch19//'.NOMA'
    innugd = ssch19//'.NUGD'
    ouvale = sdeval//'.VALE'
    oupnbn = sdeval//'.PNBN'
    oupnb2 = sdeval//'.PNB2'
    oupnco = sdeval//'.PNCO'
    oupnsp = sdeval//'.PNSP'
    oupadr = sdeval//'.PADR'
    oupcmp = sdeval//'.PCMP'
    ounoma = sdeval//'.NOMA'
    ounugd = sdeval//'.NUGD'
    ouerre = sdeval//'.ERRE'
    nabsc = sdlieu//'.ABSC'
    nrefe = sdlieu//'.REFE'
    ndesc = sdlieu//'.DESC'
    call jelira(invale, 'DOCU', cval=docu)
    call jeveuo(nrefe, 'L', arefe)
    call jeveuo(ndesc, 'L', adesc)
    call jelira(jexnum(nabsc, 1), 'LONMAX', nbnpst)
    nbmpst = nbnpst
    call jelira(inpcmp, 'LONMAX', nbtcmp)
    call jeveuo(inpcmp, 'L', aipcmp)
    call wkvect(oupcmp, 'V V I', nbtcmp, aopcmp)
!
    nbcmp = 0
    do i = 1, nbtcmp, 1
        nbcmp = nbcmp+min(zi(aipcmp+i-1), 1)
        zi(aopcmp+i-1) = zi(aipcmp+i-1)
    end do
    call wkvect(ounoma, 'V V K8', 1, adrou)
    call jeveuo(innoma, 'L', adrin)
    mailla = zk8(adrin)
    zk8(adrou) = mailla
    call wkvect(ounugd, 'V V I', 1, adrou)
    call jeveuo(innugd, 'L', ainugd)
    zi(adrou) = zi(ainugd)
    nbocer = nbmpst
!
    call jecrec(ouerre, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbocer)
    do iocer = 1, nbocer, 1
        call jecroc(jexnum(ouerre, iocer))
        call jeecra(jexnum(ouerre, iocer), 'LONMAX', nbcmp)
        call jeveuo(jexnum(ouerre, iocer), 'E', aoerre)
        do i = 1, nbcmp, 1
            zi(aoerre+i-1) = 0
        end do
    end do
!
    call jeexin(mailla//'.NOMNOE', ier)
    lnomnoe = .false.
    if (ier .ne. 0) then
        lnomnoe = .true.
    end if
    call wkvect(nnumnd, 'V V I', nbnpst, anumnd)
    do i = 1, nbnpst, 1
        if (lnomnoe) then
            call jenonu(jexnom(mailla//'.NOMNOE', zk8(adesc+i-1)), zi(anumnd+i-1))
        else
            zi(anumnd+i-1) = char8_to_int(zk8(adesc+i-1))
        end if
    end do
    call wkvect(oupadr, 'V V I', nbnpst, aopadr)
    call jeveuo(invale, 'L', aivale)
    call jeveuo(inpadr, 'L', aipadr)
!
    if (docu .eq. 'CHNO') then
        call wkvect(ouvale, 'V V R', nbcmp*nbnpst, aovale)
        zi(aopadr+1-1) = 1
        do i = 1, nbnpst-1, 1
            zi(aopadr+i+1-1) = zi(aopadr+i-1)+nbcmp
        end do
        do i = 1, nbnpst, 1
            adrin = zi(aipadr+zi(anumnd+i-1)-1)
            adrou = zi(aopadr+i-1)
            do j = 1, nbcmp, 1
                zr(aovale+adrou+j-2) = zr(aivale+adrin+j-2)
            end do
        end do
!
    else if (docu .eq. 'CHLM') then
        nindir = '&&RVECHN.TABLE.INDIR'
        call tremno(zk8(adesc+1-1), ssch19, sdemno)
        call wkvect(oupnbn, 'V V I', nbnpst, aopnbn)
        call wkvect(oupnb2, 'V V I', nbnpst, aopnb2)
        call wkvect(oupnco, 'V V I', nbnpst, aopnco)
        call wkvect(oupnsp, 'V V I', nbnpst, aopnsp)
        call wkvect(nindir, 'V V I', nbnpst, aindir)
        call jeveuo(sdemno//'.NUND', 'L', vi=nund)
        call jelira(sdemno//'.NUND', 'LONMAX', lpt)
        call jeveuo(inpnco, 'L', aipnco)
        call jeveuo(inpnsp, 'L', aipnsp)
        nnmail = sdeval//'.MAIL'
        nrepma = mailla//'.NOMMAI'
        call jeexin(nrepma, ier)
        lnommai = .false.
        if (ier .ne. 0) then
            lnommai = .true.
        end if
        do i = 1, nbnpst, 1
            pt = 1
            trouve = .false.
            n = zi(anumnd+i-1)
210         continue
            if ((.not. trouve) .and. (pt .le. lpt)) then
                if (nund(pt) .eq. n) then
                    trouve = .true.
                    call jelira(jexnum(sdemno//'.NUMA', pt), 'LONMAX', nbm)
                    zi(aindir+i-1) = pt
                    zi(aopnb2+i-1) = nbm
                    zi(aopnbn+i-1) = 1
                end if
                pt = pt+1
                goto 210
            end if
            if (.not. trouve) then
                vali = n
                valk = zk8(adesc+i-1)
                call utmess('F', 'POSTRELE_40', sk=valk, si=vali)
            end if
            call jeveuo(jexnum(sdemno//'.NUMA', pt-1), 'L', anuma)
            nsp = zi(aipnsp+zi(anuma)-1)
            nco = zi(aipnco+zi(anuma)-1)
            do j = 2, nbm, 1
                nsp = min(nsp, zi(aipnsp+zi(anuma+j-1)-1))
                nco = min(nco, zi(aipnco+zi(anuma+j-1)-1))
            end do
            zi(aopnsp+i-1) = nsp
            zi(aopnco+i-1) = nco
        end do
!
        zi(aopadr+1-1) = 1
        do i = 1, nbnpst-1, 1
            zi(aopadr+i+1-1) = zi(aopadr+i-1)+nbcmp*zi(aopnbn+i-1)*zi(aopnco+i-1)*zi(aopnsp&
                                 &+i-1)
        end do
        ilong = zi(aopadr+nbnpst-1)+nbcmp*zi(aopnbn+nbnpst-1)*zi(aopnco+nbnpst-1)*zi(aopnsp+nb&
                &npst-1)-1
        call wkvect(ouvale, 'V V R', ilong, aovale)
        call jecrec(nnmail, 'V V K8', 'NU', 'DISPERSE', 'VARIABLE', &
                    nbnpst)
        do i = 1, nbnpst, 1
            l = zi(aopnb2+i-1)
            call jecroc(jexnum(nnmail, i))
            call jeecra(jexnum(nnmail, i), 'LONMAX', l)
            call jeveuo(jexnum(nnmail, i), 'E', anmail)
            call jeveuo(jexnum(sdemno//'.NUMA', zi(aindir+i-1)), 'L', anuma)
            do j = 1, l, 1
                if (lnommai) then
                    call jenuno(jexnum(nrepma, zi(anuma+j-1)), zk8(anmail+j-1))
                else
                    zk8(anmail+j-1) = int_to_char8(zi(anuma+j-1))
                end if
            end do
        end do
        call jedetr(sdemno//'.VACP')
        call jedetr(sdemno//'.NUMA')
        call jedetr(sdemno//'.NOCP')
        call jedetr(sdemno//'.NUCP')
        call jedetr(sdemno//'.NUND')
        call jeveuo(jexnum('&CATA.GD.NOMCMP', zi(ainugd)), 'L', acmpgd)
        do i = 1, nbtcmp, 1
!
            pt = zi(aopcmp+i-1)
            if (pt .gt. 0) then
                call tremno(zk8(acmpgd+i-1), ssch19, sdemno)
                do j = 1, nbnpst, 1
                    call jeveuo(jexnum(sdemno//'.VACP', zi(aindir+j-1)), 'L', sdvacp)
                    call jeveuo(jexnum(sdemno//'.NUMA', zi(aindir+j-1)), 'L', anuma)
                    nsp = zi(aopnsp+j-1)
                    nco = zi(aopnco+j-1)
                    nbm = zi(aopnb2+j-1)
                    lnc = zi(aopadr+j-1)
                    lmc = nbcmp*nsp
!*+*                  LCC = LMC*NBM
                    lcc = lmc
                    lms = 0
                    do m = 1, nbm, 1
                        numm = zi(anuma+m-1)
                        nspm = zi(aipnsp+numm-1)
                        ncom = zi(aipnco+numm-1)
                        do k = 1, nco, 1
                            lsc = (k-1)*lcc
                            do l = 1, nsp, 1
                                indi1 = lnc-1+lsc+(l-1)*nbcmp+pt-1
                                indi2 = lms+(k-1)*nspm+l-1
                                if (zr(sdvacp+indi2) .eq. r8vide()) goto 224
                                zr(aovale+indi1) = zr(aovale+indi1)+zr(sdvacp+indi2)
224                             continue
                            end do
!
                        end do
                        lms = lms+nspm*ncom
                    end do
!
                    if (nbm .gt. 1) then
                        do k = 1, nco, 1
                            lsc = (k-1)*lcc
                            do l = 1, nsp, 1
                                indi1 = lnc-1+lsc+(l-1)*nbcmp+pt-1
                                zr(aovale+indi1) = zr(aovale+indi1)/nbm
                            end do
                        end do
                    end if
!
                end do
!
                call jedetr(sdemno//'.VACP')
                call jedetr(sdemno//'.NUMA')
                call jedetr(sdemno//'.NOCP')
                call jedetr(sdemno//'.NUCP')
                call jedetr(sdemno//'.NUND')
            end if
        end do
        call jedetr(nindir)
    else
    end if
    call jeecra(ouvale, 'DOCU', cval=docu)
    call jedetr(nnumnd)
    call jedetr(oupnb2)
    call jedema()
end subroutine
