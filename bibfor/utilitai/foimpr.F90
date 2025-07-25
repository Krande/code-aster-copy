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
subroutine foimpr(nomf, impr, iul, ind, fonins)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/foec1c.h"
#include "asterfort/foec1f.h"
#include "asterfort/foec1n.h"
#include "asterfort/foec2c.h"
#include "asterfort/foec2f.h"
#include "asterfort/foec2n.h"
#include "asterfort/fointc.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: nomf, fonins
    integer(kind=8) :: impr, iul, ind
!     ROUTINE D'IMPRESSION D'UNE FONCTION SUR UN FICHIER
!     ----------------------------------------------------------------
!
    character(len=19) :: nomfon, nomf1, listr
    character(len=24) :: prol, vale, para
    character(len=24) :: nompar, nomres, titr
    integer(kind=8) :: nbpu
    character(len=8) :: nompu
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ideb, ifin, ii, iret, ival, jval
    integer(kind=8) :: lfon, lfon1, lprol, lprol1, ltitr, lval
    integer(kind=8) :: lval1, nbfonc, nbnova, nbtitr, nbv, nbv2
    integer(kind=8) :: nbval
    real(kind=8) :: resuim, resure
    character(len=24), pointer :: nova(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    if (impr .le. 0) goto 999
    if (iul .le. 0) then
        call utmess('A', 'UTILITAI2_7')
        goto 999
    end if
    listr = fonins
    nomf1 = '&&FOIMPR'
!
!     --- NOM DE LA FONCTION A EDITER ---
    nomfon = nomf
    prol = nomfon//'.PROL'
    vale = nomfon//'.VALE'
    para = nomfon//'.PARA'
    titr = nomfon//'.TITR'
!
!     --- IMPRESSION DU TITRE ---
    write (iul, '(/,80(''-''))')
    call jeexin(titr, iret)
    if (iret .ne. 0) then
        call jeveuo(titr, 'L', ltitr)
        call jelira(titr, 'LONMAX', nbtitr)
        do i = 1, nbtitr
            write (iul, *) zk80(ltitr+i-1)
        end do
    end if
!
!     --- CAS D'UNE FONCTION "FORMULE" ---
    call jeexin(nomfon//'.NOVA', iret)
    if (iret .ne. 0 .and. ind .ne. 0) then
        call jeveuo(nomfon//'.NOVA', 'L', vk24=nova)
        call jelira(nomfon//'.NOVA', 'LONUTI', nbnova)
        if (nbnova .ne. 1) then
            call utmess('A', 'UTILITAI2_8')
            goto 999
        end if
        call jeveuo(listr//'.VALE', 'L', jval)
        call jelira(listr//'.VALE', 'LONUTI', nbval)
        nbv = 2*nbval
        call wkvect(nomf1//'.VALE', 'V V R8', nbv, lval1)
        lfon1 = lval1+nbval
        do ival = 0, nbval-1
            zr(lval1+ival) = zr(jval+ival)
            call fointe('F ', nomfon, nbnova, nova, zr(lval1+ival), &
                        zr(lfon1+ival), iret)
        end do
!
        ASSERT(lxlgut(nomf1) .le. 24)
        call wkvect(nomf1//'.PROL', 'V V K24', 6, lprol1)
        zk24(lprol1) = 'FONCTION'
        zk24(lprol1+1) = 'LIN LIN '
        zk24(lprol1+2) = nova(1)
        zk24(lprol1+3) = 'TOUTRESU'
        zk24(lprol1+4) = 'EE'
        zk24(lprol1+5) = nomf1
!
        call foec1f(iul, nomfon, zk24(lprol1), nbval, 'RIEN')
        if (impr .ge. 2) then
            ideb = 1
            ifin = min(10, nbval)
            if (impr .ge. 3) ifin = nbval
            nompar = zk24(lprol1+2)
            nomres = zk24(lprol1+3)
            call foec2f(iul, zr(lval1), nbval, ideb, ifin, &
                        nompar, nomres)
        end if
        call jedetr(nomf1//'.PROL')
        call jedetr(nomf1//'.VALE')
        goto 999
    end if
!
!     --- INFORMATIONS COMPLEMENTAIRES POUR L'EDITION ---
    call jeveuo(prol, 'L', lprol)
    nompar = zk24(lprol+2)
    nomres = zk24(lprol+3)
!
    if (zk24(lprol) .eq. 'CONSTANT' .or. zk24(lprol) .eq. 'FONCTION') then
!
!        --- NOMBRE DE VALEURS DE LA FONCTION ---
        if (ind .ne. 0) then
            call jelira(listr//'.VALE', 'LONUTI', nbval)
        else
            call jelira(vale, 'LONUTI', nbval)
            nbval = nbval/2
        end if
!
        call foec1f(iul, nomfon, zk24(lprol), nbval, 'RIEN')
        if (impr .ge. 2) then
            call jeveuo(vale, 'L', lval)
            if (ind .ne. 0) then
                call jeveuo(listr//'.VALE', 'L', jval)
                nbv2 = 2*nbval
                call wkvect(nomf1//'.VALE', 'V V R8', nbv2, lval)
                lfon = lval+nbval
                do ival = 0, nbval-1
                    zr(lval+ival) = zr(jval+ival)
                    call fointe('F ', nomfon, 1, nompar, zr(lval+ival), &
                                zr(lfon+ival), iret)
                end do
            end if
            ideb = 1
            ifin = min(10, nbval)
            if (impr .ge. 3) ifin = nbval
            call foec2f(iul, zr(lval), nbval, ideb, ifin, &
                        nompar, nomres)
            if (ind .ne. 0) then
                call jedetr(nomf1//'.PROL')
                call jedetr(nomf1//'.VALE')
            end if
        end if
!
    else if (zk24(lprol) .eq. 'NAPPE   ') then
!
        para = nomfon//'.PARA'
        call jelira(para, 'LONMAX', nbfonc)
        call foec1n(iul, nomfon, zk24(lprol), nbfonc, 'RIEN')
        if (impr .ge. 2) then
            call jeveuo(para, 'L', lval)
            ASSERT(ind .eq. 0)
            call foec2n(iul, zk24(lprol), zr(lval), vale, nbfonc, &
                        impr)
        end if
!
    else if (zk24(lprol) .eq. 'FONCT_C ') then
!
        nbpu = 1
        nompu = ' '
        call jelira(vale, 'LONUTI', nbval)
        nbval = nbval/3
        call foec1c(iul, nomfon, zk24(lprol), nbval, 'RIEN')
        if (impr .ge. 2) then
            call jeveuo(vale, 'L', lval)
            if (ind .ne. 0) then
                call jeveuo(listr//'.VALE', 'L', jval)
                call jelira(listr//'.VALE', 'LONUTI', nbval)
                nbv2 = 3*nbval
                call wkvect(nomf1//'.VALE', 'V V R8', nbv2, lval)
                lfon = lval+nbval
                ii = 0
                do ival = 0, nbval-1
                    zr(lval+ival) = zr(jval+ival)
                    call fointc('F', nomfon, nbpu, nompu, zr(lval+ival), &
                                resure, resuim, iret)
                    zr(lfon+ii) = resure
                    ii = ii+1
                    zr(lfon+ii) = resuim
                    ii = ii+1
                end do
            end if
            ideb = 1
            ifin = min(10, nbval)
            if (impr .ge. 3) ifin = nbval
            call foec2c(iul, zr(lval), nbval, ideb, ifin, &
                        nompar, nomres)
            if (ind .ne. 0) then
                call jedetr(nomf1//'.PROL')
                call jedetr(nomf1//'.VALE')
            end if
        end if
!
    else if (zk24(lprol) .eq. 'INTERPRE') then
        call utmess('A', 'UTILITAI2_10', sk=zk24(lprol))
!
    else
        call utmess('A', 'UTILITAI2_11', sk=zk24(lprol))
!
    end if
999 continue
    call jedema()
end subroutine
