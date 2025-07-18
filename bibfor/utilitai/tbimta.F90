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
subroutine tbimta(table, ifr, nparim, lipaim, formar)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: ifr, nparim
    character(len=*) :: table, formar, lipaim(*)
! ----------------------------------------------------------------------
!      IMPRESSION DE LA TABLE AU FORMAT "TABLEAU"
! ----------------------------------------------------------------------
! IN  : TABLE  : NOM D'UNE STRUCTURE "TABLE"
! IN  : IFR    : NUMERO D'UNITE LOGIQUE D'IMPRESSION
! IN  : NPARIM : NOMBRE DE PARAMETRES D'IMPRESSION
! IN  : LIPAIM : LISTE DES PARAMETRES D'IMPRESSION
! IN  : FORMAR : FORMAT D'IMPRESSION DES REELS
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, ipar, jvale, jlogq, ideb, ifin
    integer(kind=8) :: itc, itl, lgt, jcol, jlig, ifinc, il, ic
    integer(kind=8) :: lgl, vi(2), vali, iret, i1, i2, i3, i4, itc1, itc2
    integer(kind=8) :: ilon, id, if, ir, nbpara, npara, icf
    integer(kind=8) :: nblign
    real(kind=8) :: vr(2), valr, prec(2)
    complex(kind=8) :: vc(2), valc
    aster_logical :: erreur
    character(len=3) :: typec, typel, ctype
    character(len=4) :: kfin, chfin
    character(len=8) :: crit(2), form1
    character(len=16) :: inpar, knpar, formr, formd
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, lipacr(2)
    character(len=24) :: valkk(4)
    character(len=80) :: vk(2), valk
    character(len=2000) :: chaine, chainc
    integer(kind=8), pointer :: nom_para(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    nomtab = table
    crit(1) = 'EGALITE'
    crit(2) = 'EGALITE'
    prec(1) = 1.d-3
    prec(2) = 1.d-3
!
    ilon = lxlgut(formar)
    formr = '('//formar(1:ilon)//')'
    id = 0
    if = 0
    do i = 1, ilon-1
        if (formar(i:i) .eq. 'D' .or. formar(i:i) .eq. 'E' .or. formar(i:i) .eq. 'F' .or. &
            formar(i:i) .eq. 'G') then
            id = i+1
        else if (formar(i:i) .eq. '.') then
            if = i-1
        end if
    end do
    if (id .eq. if .and. id .ne. 0) then
        read (formar(id:if), '(I1)') ir
    else if (id+1 .eq. if) then
        read (formar(id:if), '(I2)') ir
    else
        ir = 12
    end if
!
    call jeveuo(nomtab//'.TBNP', 'L', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
!
!     --- ON STOCKE LES POINTEURS POUR NE PLUS FAIRE DE JEVEUO ---
!
    AS_ALLOCATE(vi=nom_para, size=nparim)
    erreur = .false.
    npara = 0
    do i = 1, nparim
        inpar = lipaim(i)
        do j = 1, nbpara
            knpar = tblp(1+4*(j-1))
            if (inpar .eq. knpar) then
                npara = npara+1
                nom_para(npara) = j
                goto 10
            end if
        end do
        erreur = .true.
        valkk(1) = inpar
        call utmess('A', 'UTILITAI6_89', sk=valkk(1))
10      continue
    end do
    if (erreur) then
        call utmess('F', 'PREPOST_60')
    end if
    if (npara .ne. 3) then
        call utmess('F', 'UTILITAI4_86')
    end if
!
    ipar = nom_para(3)
    typec = tblp(1+4*(ipar-1)+1)
    if (typec(1:1) .eq. 'I') then
        itc1 = 12
    else if (typec(1:1) .eq. 'R') then
        itc1 = ir
    else if (typec(1:1) .eq. 'C') then
        itc1 = 1+(2*ir)
    else if (typec(1:3) .eq. 'K80') then
        itc1 = 80
    else if (typec(1:3) .eq. 'K32') then
        itc1 = 32
    else if (typec(1:3) .eq. 'K24') then
        itc1 = 24
    else if (typec(1:3) .eq. 'K16') then
        itc1 = 16
    else if (typec(1:2) .eq. 'K8') then
        itc1 = 8
    end if
!
    ipar = nom_para(1)
    lipacr(2) = tblp(1+4*(ipar-1))
    typec = tblp(1+4*(ipar-1)+1)
    nomjv = tblp(1+4*(ipar-1)+2)
    nomjvl = tblp(1+4*(ipar-1)+3)
    call jeveuo(nomjv, 'L', jvale)
    call jeveuo(nomjvl, 'L', jlogq)
    if (typec(1:1) .eq. 'I') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V I', nblign, jcol)
        itc2 = 12
    else if (typec(1:1) .eq. 'R') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V R', nblign, jcol)
        itc2 = ir
    else if (typec(1:1) .eq. 'C') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V C', nblign, jcol)
        itc2 = 1+(2*ir)
    else if (typec(1:3) .eq. 'K80') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V K80', nblign, jcol)
        itc2 = 80
    else if (typec(1:3) .eq. 'K32') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V K32', nblign, jcol)
        itc2 = 32
    else if (typec(1:3) .eq. 'K24') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V K24', nblign, jcol)
        itc2 = 24
    else if (typec(1:3) .eq. 'K16') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V K16', nblign, jcol)
        itc2 = 16
    else if (typec(1:2) .eq. 'K8') then
        call wkvect('&&TBIMTA.VALE_COL', 'V V K8', nblign, jcol)
        itc2 = 8
    end if
    itc = max(itc1, itc2)
    ic = 0
    do i = 1, nblign
        if (zi(jlogq+i-1) .eq. 1) then
            if (typec(1:1) .eq. 'I') then
                do j = 1, ic
                    if (zi(jcol+j-1) .eq. zi(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zi(jcol+ic-1) = zi(jvale+i-1)
            else if (typec(1:1) .eq. 'R') then
                do j = 1, ic
                    if (zr(jcol+j-1) .eq. zr(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zr(jcol+ic-1) = zr(jvale+i-1)
            else if (typec(1:1) .eq. 'C') then
                do j = 1, ic
                    if (zc(jcol+j-1) .eq. zc(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zc(jcol+ic-1) = zc(jvale+i-1)
            else if (typec(1:3) .eq. 'K80') then
                do j = 1, ic
                    if (zk80(jcol+j-1) .eq. zk80(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zk80(jcol+ic-1) = zk80(jvale+i-1)
            else if (typec(1:3) .eq. 'K32') then
                do j = 1, ic
                    if (zk32(jcol+j-1) .eq. zk32(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zk32(jcol+ic-1) = zk32(jvale+i-1)
            else if (typec(1:3) .eq. 'K24') then
                do j = 1, ic
                    if (zk24(jcol+j-1) .eq. zk24(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zk24(jcol+ic-1) = zk24(jvale+i-1)
            else if (typec(1:3) .eq. 'K16') then
                do j = 1, ic
                    if (zk16(jcol+j-1) .eq. zk16(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zk16(jcol+ic-1) = zk16(jvale+i-1)
            else if (typec(1:2) .eq. 'K8') then
                do j = 1, ic
                    if (zk8(jcol+j-1) .eq. zk8(jvale+i-1)) goto 100
                end do
                ic = ic+1
                zk8(jcol+ic-1) = zk8(jvale+i-1)
            end if
        end if
100     continue
    end do
!
    ipar = nom_para(2)
    lipacr(1) = tblp(1+4*(ipar-1))
    typel = tblp(1+4*(ipar-1)+1)
    nomjv = tblp(1+4*(ipar-1)+2)
    nomjvl = tblp(1+4*(ipar-1)+3)
    lgl = lxlgut(lipacr(1))
    call jeveuo(nomjv, 'L', jvale)
    call jeveuo(nomjvl, 'L', jlogq)
    if (typel(1:1) .eq. 'I') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V I', nblign, jlig)
        itl = 12
    else if (typel(1:1) .eq. 'R') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V R', nblign, jlig)
        itl = ir
    else if (typel(1:1) .eq. 'C') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V C', nblign, jlig)
        itl = 1+(2*ir)
    else if (typel(1:3) .eq. 'K80') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V K80', nblign, jlig)
        itl = 80
    else if (typel(1:3) .eq. 'K32') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V K32', nblign, jlig)
        itl = 32
    else if (typel(1:3) .eq. 'K24') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V K24', nblign, jlig)
        itl = 24
    else if (typel(1:3) .eq. 'K16') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V K16', nblign, jlig)
        itl = 16
    else if (typel(1:2) .eq. 'K8') then
        call wkvect('&&TBIMTA.VALE_LIG', 'V V K8', nblign, jlig)
        itl = 8
    end if
    il = 0
    do i = 1, nblign
        if (zi(jlogq+i-1) .eq. 1) then
            if (typel(1:1) .eq. 'I') then
                do j = 1, il
                    if (zi(jlig+j-1) .eq. zi(jvale+i-1)) goto 200
                end do
                il = il+1
                zi(jlig+il-1) = zi(jvale+i-1)
            else if (typel(1:1) .eq. 'R') then
                do j = 1, il
                    if (zr(jlig+j-1) .eq. zr(jvale+i-1)) goto 200
                end do
                il = il+1
                zr(jlig+il-1) = zr(jvale+i-1)
            else if (typel(1:1) .eq. 'C') then
                do j = 1, il
                    if (zc(jlig+j-1) .eq. zc(jvale+i-1)) goto 200
                end do
                il = il+1
                zc(jlig+il-1) = zc(jvale+i-1)
            else if (typel(1:3) .eq. 'K80') then
                do j = 1, il
                    if (zk80(jlig+j-1) .eq. zk80(jvale+i-1)) goto 200
                end do
                il = il+1
                zk80(jlig+il-1) = zk80(jvale+i-1)
            else if (typel(1:3) .eq. 'K32') then
                do j = 1, il
                    if (zk32(jlig+j-1) .eq. zk32(jvale+i-1)) goto 200
                end do
                il = il+1
                zk32(jlig+il-1) = zk32(jvale+i-1)
            else if (typel(1:3) .eq. 'K24') then
                do j = 1, il
                    if (zk24(jlig+j-1) .eq. zk24(jvale+i-1)) goto 200
                end do
                il = il+1
                zk24(jlig+il-1) = zk24(jvale+i-1)
            else if (typel(1:3) .eq. 'K16') then
                do j = 1, il
                    if (zk16(jlig+j-1) .eq. zk16(jvale+i-1)) goto 200
                end do
                il = il+1
                zk16(jlig+il-1) = zk16(jvale+i-1)
            else if (typel(1:2) .eq. 'K8') then
                do j = 1, il
                    if (zk8(jlig+j-1) .eq. zk8(jvale+i-1)) goto 200
                end do
                il = il+1
                zk8(jlig+il-1) = zk8(jvale+i-1)
            end if
        end if
200     continue
    end do
    lgt = lgl+itl+1
!
    icf = ic
    chainc = ' '
    ideb = lgt+1
    ifin = ideb+3
    chainc(ideb:ifin) = ' ! '
    ideb = ifin+1
    do i = 1, icf
        ifin = ideb+itc-1
        if (ifin .gt. 2000) then
            icf = i-1
            ifin = ideb
            goto 302
        end if
        if (typec(1:1) .eq. 'I') then
            write (chainc(ideb:ifin), '(I12)') zi(jcol+i-1)
        else if (typec(1:1) .eq. 'R') then
            write (chainc(ideb:ifin), formr) zr(jcol+i-1)
        else if (typec(1:1) .eq. 'C') then
            write (chainc(ideb:ifin), formr) zc(jcol+i-1)
        else if (typec(1:3) .eq. 'K80') then
            chainc(ideb:ifin) = zk80(jcol+i-1)
        else if (typec(1:3) .eq. 'K32') then
            chainc(ideb:ifin) = zk32(jcol+i-1)
        else if (typec(1:3) .eq. 'K24') then
            chainc(ideb:ifin) = zk24(jcol+i-1)
        else if (typec(1:3) .eq. 'K16') then
            chainc(ideb:ifin) = zk16(jcol+i-1)
        else if (typec(1:2) .eq. 'K8') then
            chainc(ideb:ifin) = zk8(jcol+i-1)
        end if
        ideb = ifin+2
    end do
302 continue
    ifinc = ifin+2
!
    call codent(ifinc, 'D', kfin)
    formd = '('//kfin//'(''-''))'
!
    ipar = nom_para(3)
    inpar = tblp(1+4*(ipar-1))
    nomjv = tblp(1+4*(ipar-1)+2)
    nomjvl = tblp(1+4*(ipar-1)+3)
    chaine = ' '
    chaine(1:lgt) = inpar
    ifin = lgt+3+24
    chaine(lgt+1:ifin) = ' ! '//lipacr(2)
!
    call codent(ifin, 'G', chfin)
    form1 = '(A'//chfin//')'
    write (ifr, form1) chaine(1:ifin)
!
    call codent(ifinc, 'G', chfin)
    form1 = '(A'//chfin//')'
    write (ifr, form1) chainc(1:ifinc)
!
    write (ifr, formd)
!
    do j = 1, il
        i1 = 1
        i2 = 1
        i3 = 1
        i4 = 1
        chaine = ' '
        if (j .eq. 1) chaine(1:lgl) = lipacr(1)
        ideb = lgl+2
        ifin = lgt
        if (typel(1:1) .eq. 'I') then
            write (chaine(ideb:ifin), '(I12)') zi(jlig+j-1)
            vi(i1) = zi(jlig+j-1)
            i1 = i1+1
        else if (typel(1:1) .eq. 'R') then
            write (chaine(ideb:ifin), formr) zr(jlig+j-1)
            vr(i2) = zr(jlig+j-1)
            i2 = i2+1
        else if (typel(1:1) .eq. 'C') then
            write (chaine(ideb:ifin), formr) zc(jlig+j-1)
            vr(i3) = dble(zc(jlig+j-1))
            i3 = i3+1
        else if (typel(1:3) .eq. 'K80') then
            chaine(ideb:ifin) = zk80(jlig+j-1)
            vk(i4) = zk80(jlig+j-1)
            i4 = i4+1
        else if (typel(1:3) .eq. 'K32') then
            chaine(ideb:ifin) = zk32(jlig+j-1)
            vk(i4) = zk32(jlig+j-1)
            i4 = i4+1
        else if (typel(1:3) .eq. 'K24') then
            chaine(ideb:ifin) = zk24(jlig+j-1)
            vk(i4) = zk24(jlig+j-1)
            i4 = i4+1
        else if (typel(1:3) .eq. 'K16') then
            chaine(ideb:ifin) = zk16(jlig+j-1)
            vk(i4) = zk16(jlig+j-1)
            i4 = i4+1
        else if (typel(1:2) .eq. 'K8') then
            chaine(ideb:ifin) = zk8(jlig+j-1)
            vk(i4) = zk8(jlig+j-1)
            i4 = i4+1
        end if
        ideb = ifin+1
        ifin = ideb+3
        chaine(ideb:ifin) = ' ! '
        ideb = ifin+1
!
        do k = 1, icf
            if (typec(1:1) .eq. 'I') then
                vi(i1) = zi(jcol+k-1)
            else if (typec(1:1) .eq. 'R') then
                vr(i2) = zr(jcol+k-1)
            else if (typec(1:1) .eq. 'C') then
                vc(i3) = zc(jcol+k-1)
            else if (typec(1:3) .eq. 'K80') then
                vk(i4) = zk80(jcol+k-1)
            else if (typec(1:3) .eq. 'K32') then
                vk(i4) = zk32(jcol+k-1)
            else if (typec(1:3) .eq. 'K24') then
                vk(i4) = zk24(jcol+k-1)
            else if (typec(1:3) .eq. 'K16') then
                vk(i4) = zk16(jcol+k-1)
            else if (typec(1:2) .eq. 'K8') then
                vk(i4) = zk8(jcol+k-1)
            end if
            call tbliva(table, 2, lipacr, vi, vr, &
                        vc, vk, crit, prec, inpar, &
                        ctype, vali, valr, valc, valk, &
                        iret)
            if (iret .eq. 0) then
                ifin = ideb+itc-1
                if (ctype(1:1) .eq. 'I') then
                    write (chaine(ideb:ifin), '(I12)') vali
                else if (ctype(1:1) .eq. 'R') then
                    write (chaine(ideb:ifin), formr) valr
                else if (ctype(1:1) .eq. 'C') then
                    write (chaine(ideb:ifin), formr) valc
                else if (ctype(1:1) .eq. 'K') then
                    chaine(ideb:ifin) = valk
                end if
                ideb = ifin+2
            else if (iret .eq. 2) then
                ifin = ideb+itc-1
                chaine(ideb:ifin) = '    -    '
                ideb = ifin+2
            else
                valkk(1) = inpar
                valkk(2) = lipacr(1)
                valkk(3) = lipacr(2)
                call utmess('F', 'UTILITAI6_99', nk=3, valk=valkk)
            end if
        end do
!
        call codent(ifin, 'G', chfin)
        form1 = '(A'//chfin//')'
        write (ifr, form1) chaine(1:ifin)
!
    end do
!
    write (ifr, formd)
!
    if (icf .ne. ic) then
        call utmess('A', 'UTILITAI4_84')
    end if
!
    AS_DEALLOCATE(vi=nom_para)
    call jedetr('&&TBIMTA.VALE_COL')
    call jedetr('&&TBIMTA.VALE_LIG')
!
    call jedema()
!
end subroutine
