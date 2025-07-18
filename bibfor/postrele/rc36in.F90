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
subroutine rc36in(noma, nbma, listma, chindi)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/rc36zz.h"
#include "asterfort/reliem.h"
!
    integer(kind=8) :: nbma, listma(*)
    character(len=8) :: noma
    character(len=24) :: chindi
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!     RECUPERATION DES DONNEES DE "INDI_SIGM"
!
! IN  : NOMA   : MAILLAGE
! IN  : NBMA   : NOMBRE DE MAILLES D'ANALYSE
! IN  : LISTMA : LISTE DES MAILLES D'ANALYSE
! OUT : CHINDI : CHAM_ELEM DE TYPE ELNO D'INDICES DE CONTRAINTES
!     ------------------------------------------------------------------
!
    integer(kind=8) :: n1, n2, nbindi, iocc, nbcmp, decal, ipt, icmp, iad, nbpt
    integer(kind=8) :: jconx2, in, im, ima, ino, nbnoeu, jnoeu, nbmail
    integer(kind=8) :: jmail, nbtou, im1
    parameter(nbcmp=7)
    real(kind=8) :: vale(nbcmp)
    character(len=8) :: k8b, nomgd, type
    character(len=8) :: motcls(2), typmcs(2), motcln(2), typmcn(2)
    character(len=16) :: motclf, nocmp(nbcmp)
    character(len=24) :: mesmai, mesnoe
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: cesd(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    motclf = 'INDI_SIGM'
!
    mesmai = 'RC36IN.MES_MAILLES'
    motcls(1) = 'GROUP_MA'
    motcls(2) = 'MAILLE'
    typmcs(1) = 'GROUP_MA'
    typmcs(2) = 'MAILLE'
    mesnoe = 'RC36IN.MES_NOEUDS'
    motcln(1) = 'GROUP_NO'
    motcln(2) = 'NOEUD'
    typmcn(1) = 'GROUP_NO'
    typmcn(2) = 'NOEUD'
!
    call getfac(motclf, nbindi)
!
    nomgd = 'RCCM_R'
    nocmp(1) = 'C1'
    nocmp(2) = 'C2'
    nocmp(3) = 'C3'
    nocmp(4) = 'K1'
    nocmp(5) = 'K2'
    nocmp(6) = 'K3'
    nocmp(7) = 'TYPE'
!
    call rc36zz(noma, nomgd, nbcmp, nocmp, nbma, &
                listma, chindi)
!
    call jeveuo(chindi(1:19)//'.CESD', 'L', vi=cesd)
    call jeveuo(chindi(1:19)//'.CESV', 'E', vr=cesv)
!
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
!
    do iocc = 1, nbindi, 1
!
        call getvr8(motclf, 'C1', iocc=iocc, scal=vale(1), nbret=n1)
        call getvr8(motclf, 'C2', iocc=iocc, scal=vale(2), nbret=n1)
        call getvr8(motclf, 'C3', iocc=iocc, scal=vale(3), nbret=n1)
        call getvr8(motclf, 'K1', iocc=iocc, scal=vale(4), nbret=n1)
        call getvr8(motclf, 'K2', iocc=iocc, scal=vale(5), nbret=n1)
        call getvr8(motclf, 'K3', iocc=iocc, scal=vale(6), nbret=n1)
!
        call getvtx(motclf, 'TYPE_ELEM_STANDARD', iocc=iocc, scal=type, nbret=n1)
        if (n1 .eq. 0) then
            vale(7) = 0.d0
        else
            if (type(1:3) .eq. 'DRO') vale(7) = 10.d0
            if (type(1:3) .eq. 'COU') vale(7) = 20.d0
            if (type(1:3) .eq. 'TRN') vale(7) = 30.d0
            if (type(1:3) .eq. 'TEE') vale(7) = 40.d0
        end if
!
        call getvtx(motclf, 'GROUP_NO', iocc=iocc, nbval=0, nbret=n1)
        call getvtx(motclf, 'NOEUD', iocc=iocc, nbval=0, nbret=n2)
        if (n1+n2 .ne. 0) then
            call reliem(' ', noma, 'NU_NOEUD', motclf, iocc, &
                        2, motcln, typmcn, mesnoe, nbnoeu)
            call jeveuo(mesnoe, 'L', jnoeu)
        else
            nbnoeu = 0
        end if
!
        call getvtx(motclf, 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            do im = 1, nbma
                ima = listma(im)
                nbpt = cesd(5+4*(ima-1)+1)
                decal = cesd(5+4*(ima-1)+4)
                do ipt = 1, nbpt
                    do icmp = 1, nbcmp
                        iad = decal+(ipt-1)*nbcmp+icmp
                        cesv(iad) = vale(icmp)
                    end do
                end do
            end do
!
        else
            call reliem(' ', noma, 'NU_MAILLE', motclf, iocc, &
                        2, motcls, typmcs, mesmai, nbmail)
            call jeveuo(mesmai, 'L', jmail)
!
            if (nbnoeu .eq. 0) then
                do im = 1, nbmail
                    ima = zi(jmail+im-1)
                    do im1 = 1, nbma
                        if (listma(im1) .eq. ima) goto 204
                    end do
                    goto 200
204                 continue
                    nbpt = cesd(5+4*(ima-1)+1)
                    decal = cesd(5+4*(ima-1)+4)
                    do ipt = 1, nbpt
                        ino = connex(zi(jconx2+ima-1)+ipt-1)
                        do icmp = 1, nbcmp
                            iad = decal+(ipt-1)*nbcmp+icmp
                            cesv(iad) = vale(icmp)
                        end do
                    end do
200                 continue
                end do
            else
                do im = 1, nbmail
                    ima = zi(jmail+im-1)
                    do im1 = 1, nbma
                        if (listma(im1) .eq. ima) goto 304
                    end do
                    goto 300
304                 continue
                    nbpt = cesd(5+4*(ima-1)+1)
                    decal = cesd(5+4*(ima-1)+4)
                    do ipt = 1, nbpt
                        ino = connex(zi(jconx2+ima-1)+ipt-1)
                        do in = 1, nbnoeu
                            if (zi(jnoeu+in-1) .eq. ino) then
                                do icmp = 1, nbcmp
                                    iad = decal+(ipt-1)*nbcmp+icmp
                                    cesv(iad) = vale(icmp)
                                end do
                                goto 310
                            end if
                        end do
310                     continue
                    end do
300                 continue
                end do
                call jedetr(mesnoe)
            end if
            call jedetr(mesmai)
        end if
!
    end do
!
    call jedema()
end subroutine
