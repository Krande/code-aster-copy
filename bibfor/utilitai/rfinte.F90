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
subroutine rfinte(ispec)
!     OPERATEUR "RECU_FONCTION"   MOT CLE "INTE_SPEC"
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/foattr.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: ispec
    integer(kind=8) :: nbval
    integer(kind=8) :: lpro, i1, indi
    integer(kind=8) :: i, kvale
    integer(kind=8) :: ifm, niv, n2, n3, n4, n5, mxval, numi, numj, nbfreq, ifreq, num
    integer(kind=8) :: lnumi, lnumj, lvale, lfreq, lnum
    integer(kind=8) :: lnoei, lnoej, lcmpi, lcmpj
    character(len=8) :: nospec, noei, noej, cmpi, cmpj
    character(len=16) :: nomcmd, typcon, noch, nocham, abscisse
    character(len=19) :: nomfon
    character(len=24) :: chnumi, chnumj, chfreq, chvale, chnum
    character(len=24) :: chnoei, chnoej, chcmpi, chcmpj
    character(len=24) :: paray
    aster_logical :: indice
    character(len=16), pointer :: refe(:) => null()
!
! DEB------------------------------------------------------------------
!
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomfon, typcon, nomcmd)
!
!
    call getvtx(' ', 'NOM_CHAM', scal=nocham, nbret=n4)
!
    nospec = ispec
    call jeveuo(nospec//'.REFE', 'L', vk16=refe)
    noch = refe(1)
    abscisse = refe(3)
    if (n4 .eq. 0) then
        paray = 'DSP'
    else
        if (nocham .ne. noch) then
            call utmess('F', 'UTILITAI_55', sk=nocham)
        else
            paray = noch
        end if
    end if
!
    chfreq = nospec//'.DISC'
    call jeveuo(chfreq, 'L', lfreq)
!
    call getvtx(' ', 'NOEUD_I', nbval=0, nbret=n2)
    call getvis(' ', 'NUME_ORDRE_I', nbval=0, nbret=n3)
    call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=n5)
!
    indice = .false.
    indi = 0
    if (n2 .lt. 0) then
        call getvtx(' ', 'NOEUD_I', scal=noei, nbret=n4)
        call getvtx(' ', 'NOEUD_J', nbval=0, nbret=n4)
        if (n4 .eq. 0) then
            noej = noei
        else
            call getvtx(' ', 'NOEUD_J', scal=noej, nbret=n4)
        end if
        call getvtx(' ', 'NOM_CMP_I', scal=cmpi, nbret=n4)
        call getvtx(' ', 'NOM_CMP_J', nbval=0, nbret=n4)
        if (n4 .eq. 0) then
            cmpj = cmpi
        else
            call getvtx(' ', 'NOM_CMP_J', scal=cmpj, nbret=n4)
        end if
        chnoei = nospec//'.NOEI'
        chnoej = nospec//'.NOEJ'
        chcmpi = nospec//'.CMPI'
        chcmpj = nospec//'.CMPJ'
        call jeveuo(chnoei, 'L', lnoei)
        call jeveuo(chnoej, 'L', lnoej)
        call jeveuo(chcmpi, 'L', lcmpi)
        call jeveuo(chcmpj, 'L', lcmpj)
        call jelira(chnoei, 'LONMAX', mxval)
        do i1 = 1, mxval
            if ((zk8(lnoei-1+i1) .eq. noei) .and. (zk8(lnoej-1+i1) .eq. noej) .and. &
                (zk8(lcmpi-1+i1) .eq. cmpi) .and. (zk8(lcmpj-1+i1) .eq. cmpj)) then
                indi = i1
                indice = .true.
            end if
        end do
    else if (n3 .lt. 0) then
        call getvis(' ', 'NUME_ORDRE_I', scal=numi, nbret=n4)
        call getvis(' ', 'NUME_ORDRE_J', nbval=0, nbret=n4)
        if (n4 .eq. 0) then
            numj = numi
        else
            call getvis(' ', 'NUME_ORDRE_J', scal=numj, nbret=n4)
        end if
        chnumi = nospec//'.NUMI'
        chnumj = nospec//'.NUMJ'
        call jeveuo(chnumi, 'L', lnumi)
        call jeveuo(chnumj, 'L', lnumj)
        call jelira(chnumi, 'LONMAX', mxval)
        do i1 = 1, mxval
            if ((zi(lnumi-1+i1) .eq. numi) .and. (zi(lnumj-1+i1) .eq. numj)) then
                indi = i1
                indice = .true.
            end if
        end do
    else if (n5 .lt. 0) then
        call getvis(' ', 'NUME_ORDRE', scal=num, nbret=n4)
        chnum = nospec//'.NUME_ORDRE'
        call jeveuo(chnum, 'L', lnum)
        call jelira(chnum, 'LONMAX', mxval)
        do i1 = 1, mxval
            if ((zi(lnum-1+i1) .eq. num) .and. (.not. indice)) then
                indi = i1
                indice = .true.
            end if
        end do
    end if
!
    if (.not. indice) then
        call utmess('F', 'UTILITAI4_53')
    end if
!
    chfreq = nospec//'.DISC'
    call jelira(chfreq, 'LONMAX', nbfreq)
    call jeveuo(chfreq, 'L', ifreq)
!
    chvale = nospec//'.VALE'
    call jeveuo(jexnum(chvale, indi), 'L', lvale)
    call jelira(jexnum(chvale, indi), 'LONMAX', nbval)
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
    zk24(lpro+1) = 'LIN LIN '
    zk24(lpro+2) = abscisse
    zk24(lpro+3) = paray
    zk24(lpro+4) = 'LL      '
    zk24(lpro+5) = nomfon
!
! --- SURCHARGE EVENTUELLE DU .PROL
!
    call foattr(' ', 1, nomfon)
!
    if (nbval .eq. nbfreq) then
        zk24(lpro) = 'FONCTION'
        call wkvect(nomfon//'.VALE', 'G V R', 2*nbfreq, kvale)
        do i = 1, nbfreq
            zr(kvale+i-1) = zr(ifreq+i-1)
            zr(kvale+nbfreq+i-1) = zr(lvale+i-1)
        end do
    else
        zk24(lpro) = 'FONCT_C'
        call wkvect(nomfon//'.VALE', 'G V R', 3*nbfreq, kvale)
        do i = 1, nbfreq
            zr(kvale+i-1) = zr(ifreq+i-1)
            zr(kvale+nbfreq+2*(i-1)) = zr(lvale+2*(i-1))
            zr(kvale+nbfreq+2*(i-1)+1) = zr(lvale+2*(i-1)+1)
        end do
    end if
!
    call jedema()
end subroutine
