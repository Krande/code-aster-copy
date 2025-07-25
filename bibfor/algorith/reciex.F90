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
subroutine reciex(intexc, iderex, nindex, nnoeex, ncmpex, &
                  nvasex, graexc, excmod, napexc)
!    C. DUVAL
!-----------------------------------------------------------------------
!  BUT: RECUPERER LES INFORMATIONS DE TYPE EXCITATION POUR
    implicit none
!        LE CALCUL DYNAMIQUE ALEATOIRE
!
! INTEXC   /OUT/: NOM DE L INTERSPECTRE  EXCITATION
! IDEREX   /OUT/: ORDRE DE DERIVATION
! NINDEX   /OUT/: NOMBRE  D INDICES RECUPERES
! NNOEEX   /OUT/: NOMBRE DE NOEUDS DONNES EN EXCITATION
! NCMPEX   /OUT/: NOMBRE DE CMP DONNES EN EXCITATION
! NVASEX   /OUT/: NOMBRE DE VECTEURS ASSEMBLES DONNES EN EXCITATION
! GRAEXC  /OUT/ : GRANDEUR EXCITATION
! EXCMOD  /OUT/ : TYPE MODAL
! NAPEXC  /OUT/ : NOMBRE D APPUI EXCITATION (NOEUDS OU VECTASS)
!
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, ibid1, iderex, ij2, ilcmpi
    integer(kind=8) :: ilcpex, ilfex, ilindi, illex, ilnoex
    integer(kind=8) :: ilvaex, ivite, napexc, ncmpex
    integer(kind=8) :: ndim, nindex, nnoeex, nvasex
    integer(kind=8) :: vali(2)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, iret
    character(len=4) :: excmod
    character(len=8) :: intexc
    character(len=16) :: graexc
    character(len=24) :: chnumi, chnumj, chnoei, chnoej, chcmpi, chcmpj, chvale
    character(len=24) :: chfreq
    character(len=24) :: valk(5)
!
    aster_logical :: lindi, exiind
    integer(kind=8) :: lnumi, lnumj, mxval, num, lcmpi, lcmpj
    integer(kind=8) :: nbfreq, ifreq
!
    call getvid('EXCIT', 'INTE_SPEC', iocc=1, scal=intexc, nbret=ibid)
!
    call getvis('EXCIT', 'DERIVATION', iocc=1, scal=iderex, nbret=ibid)
!
    call getvis('EXCIT', 'NUME_ORDRE_I', iocc=1, nbval=0, nbret=nindex)
    if (nindex .ne. 0) then
        lindi = .true.
        nindex = -nindex
        call wkvect('&&RECIEX.INDI_I', 'V V I', nindex, ilindi)
        call getvis('EXCIT', 'NUME_ORDRE_I', iocc=1, nbval=nindex, vect=zi(ilindi), &
                    nbret=ibid)
    else
        call getvtx('EXCIT', 'NOEUD_I', iocc=1, nbval=0, nbret=nindex)
        lindi = .false.
        nindex = -nindex
        call wkvect('&&RECIEX.INDI_I', 'V V K8', nindex, ilindi)
        call wkvect('&&RECIEX.CMP_I', 'V V K8', nindex, ilcmpi)
        call getvtx('EXCIT', 'NOEUD_I', iocc=1, nbval=nindex, vect=zk8(ilindi), &
                    nbret=ibid)
        call getvtx('EXCIT', 'NOM_CMP_I', iocc=1, nbval=nindex, vect=zk8(ilcmpi), &
                    nbret=ibid)
    end if
    call getvis('EXCIT', 'NUME_VITE_FLUI', iocc=1, scal=ivite, nbret=ibid)
!
    ndim = nindex*(nindex+1)/2
    call wkvect('&&OP0131.LIADRFEX1', 'V V I', ndim, ilfex)
    call wkvect('&&OP0131.LIADRLEX1', 'V V I', ndim+2, illex)
!
    chfreq = intexc//'.DISC'
    call jelira(chfreq, 'LONMAX', nbfreq)
    call jeveuo(chfreq, 'L', ifreq)
    zi(illex) = nbfreq
    zi(illex+ndim+1) = ifreq
    chvale = intexc//'.VALE'
!
!     VERIFICATIONS EXISTENCE PARAMETRES DE LA SD
    if (lindi) then
        chnumi = intexc//'.NUMI'
        chnumj = intexc//'.NUMJ'
        call jeveuo(chnumi, 'L', lnumi)
        call jeveuo(chnumj, 'L', lnumj)
        call jelira(chnumi, 'LONMAX', mxval)
        do i1 = 1, nindex
            do i2 = i1, nindex
                ij2 = (i2*(i2-1))/2+i1
                exiind = .false.
                do num = 1, mxval
                    if (( &
                        (zi(lnumi-1+num) .eq. zi(ilindi-1+i1)) .and. &
                        (zi(lnumj-1+num) .eq. zi(ilindi-1+i2)) &
                        ) &
                        .or. &
                        ( &
                        (zi(lnumi-1+num) .eq. zi(ilindi-1+i2)) .and. &
                        (zi(lnumj-1+num) .eq. zi(ilindi-1+i1)) &
                        )) then
                        exiind = .true.
                        call jeveuo(jexnum(chvale, num), 'L', zi(ilfex-1+ij2))
                        call jelira(jexnum(chvale, num), 'LONMAX', zi(illex+ij2))
                    end if
                end do
                if (.not. exiind) then
                    valk(1) = intexc
                    vali(1) = zi(ilindi-1+i1)
                    vali(2) = zi(ilindi-1+i2)
                    call utmess('F', 'PREPOST3_84', sk=valk(1), ni=2, vali=vali)
                end if
            end do
        end do
    else
        chnoei = intexc//'.NOEI'
        chnoej = intexc//'.NOEJ'
        chcmpi = intexc//'.CMPI'
        chcmpj = intexc//'.CMPJ'
        call jeveuo(chnoei, 'L', lnumi)
        call jeveuo(chnoej, 'L', lnumj)
        call jeveuo(chcmpi, 'L', lcmpi)
        call jeveuo(chcmpj, 'L', lcmpj)
        call jelira(chnoei, 'LONMAX', mxval)
        do i1 = 1, nindex
            do i2 = i1, nindex
                ij2 = (i2*(i2-1))/2+i1
                exiind = .false.
                do num = 1, mxval
                    if (( &
                        (zk8(lnumi-1+num) .eq. zk8(ilindi-1+i1)) .and. &
                        (zk8(lnumj-1+num) .eq. zk8(ilindi-1+i2)) .and. &
                        (zk8(lcmpi-1+num) .eq. zk8(ilcmpi-1+i1)) .and. &
                        (zk8(lcmpj-1+num) .eq. zk8(ilcmpi-1+i2)) &
                        ) &
                        .or. &
                        ( &
                        (zk8(lnumi-1+num) .eq. zk8(ilindi-1+i2)) .and. &
                        (zk8(lnumj-1+num) .eq. zk8(ilindi-1+i1)) .and. &
                        (zk8(lcmpi-1+num) .eq. zk8(ilcmpi-1+i2)) .and. &
                        (zk8(lcmpj-1+num) .eq. zk8(ilcmpi-1+i1)) &
                        )) then
                        exiind = .true.
                        call jeveuo(jexnum(chvale, num), 'L', zi(ilfex-1+ij2))
                        call jelira(jexnum(chvale, num), 'LONMAX', zi(illex+ij2))
                    end if
                end do
                if (.not. exiind) then
                    valk(1) = zk8(ilindi-1+i1)
                    valk(2) = zk8(ilcmpi-1+i1)
                    valk(3) = zk8(ilindi-1+i2)
                    valk(4) = zk8(ilcmpi-1+i2)
                    valk(5) = intexc
                    call utmess('F', 'PREPOST3_85', nk=5, valk=valk)
                end if
            end do
        end do
    end if
!
!----TYPE MODAL ('NON' PAR DEFAUT)
!
    call getvtx('EXCIT', 'MODAL', iocc=1, scal=excmod, nbret=ibid)
    if (excmod .eq. 'OUI') napexc = nindex
!
!----GRANDEUR   (DEPL_R PAR DEFAUT)
!
    call getvtx('EXCIT', 'GRANDEUR', iocc=1, scal=graexc, nbret=ibid)
!
!---NOEUDS APPUIS
!
    call getvtx('EXCIT', 'NOEUD', iocc=1, nbval=0, nbret=nnoeex)
    nnoeex = -nnoeex
    if (nnoeex .ne. 0) then
        napexc = nnoeex
        call wkvect('&&OP0131.LISTENOEEXC', 'V V K8', nnoeex, ilnoex)
        call getvtx('EXCIT', 'NOEUD', iocc=1, nbval=nnoeex, vect=zk8(ilnoex), &
                    nbret=ibid)
    end if
!
!---CMP APPUIS
!
    call getvtx('EXCIT', 'NOM_CMP', iocc=1, nbval=0, nbret=ncmpex)
    ncmpex = -ncmpex
    if (ncmpex .ne. 0) then
        call wkvect('&&OP0131.LISTECMPEXC', 'V V K8', ncmpex, ilcpex)
        call getvtx('EXCIT', 'NOM_CMP', iocc=1, nbval=ncmpex, vect=zk8(ilcpex), &
                    nbret=ibid)
    end if
!
!---VECTEURS ASSEMBLES
!
    call getvid('EXCIT', 'CHAM_NO', iocc=1, nbval=0, nbret=nvasex)
    nvasex = -nvasex
    if (nvasex .ne. 0) then
        napexc = nvasex
        graexc = 'EFFO'
        call wkvect('&&OP0131.LVECTASSEXC', 'V V K8', nvasex, ilvaex)
        call getvid('EXCIT', 'CHAM_NO', iocc=1, nbval=nvasex, vect=zk8(ilvaex), &
                    nbret=ibid1)
    end if
!
    if (graexc .eq. 'EFFO') iderex = 0
!
    call jedetr('&&RECIEX.INDI_I')
    call jeexin('&&RECIEX.CMP_I', iret)
    if (iret .ne. 0) then
        call jedetr('&&RECIEX.CMP_I')
    end if
!
end subroutine
