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
subroutine tbimex(table, ifr, nparim, lipaim, formaz, &
                  formar)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: ifr, nparim
    character(len=*) :: table, lipaim(*), formaz, formar
!      IMPRESSION D'UNE TABLE AU FORMAT "EXCEL" OU "AGRAF"
! ----------------------------------------------------------------------
! IN  : TABLE  : NOM D'UNE STRUCTURE "TABLE"
! IN  : IFR    : NUMERO D'UNITE LOGIQUE D'IMPRESSION
! IN  : NPARIM : NOMBRE DE PARAMETRES D'IMPRESSION
! IN  : LIPAIM : LISTE DES PARAMETRES D'IMPRESSION
! IN  : FORMAT : FORMAT D'IMPRESSION DE LA TABLE
! IN  : FORMAR : FORMAT D'IMPRESSION DES REELS
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, ipar, jvale, jlogq, ideb, ifin
    integer(kind=8) :: nblign, ilon, ilm, id, if, ir, ilmp, iaj, nbpara, npara
    integer(kind=8) :: nparaf
    aster_logical :: erreur
    character(len=1) :: bacs
    character(len=3) :: type
    character(len=4) :: chfin
    character(len=8) :: format, form1
    character(len=16) :: formr
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, knpar
    character(len=24) :: valk
    character(len=2000) :: chaine, chain2
    integer(kind=8), pointer :: log_para(:) => null()
    integer(kind=8), pointer :: nom_para(:) => null()
    integer(kind=8), pointer :: val_para(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    nomtab = table
    format = formaz
    bacs = char(92)
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
!     --- ON RECHERCHE LA LONGUEUR LA PLUS LONGUE ---
!     --- ON STOCKE LES POINTEURS POUR NE PLUS FAIRE DE JEVEUO ---
!
    AS_ALLOCATE(vi=nom_para, size=nparim)
    AS_ALLOCATE(vi=val_para, size=nparim)
    AS_ALLOCATE(vi=log_para, size=nparim)
    erreur = .false.
    npara = 0
    ilmp = 0
    do i = 1, nparim
        inpar = lipaim(i)
        do j = 1, nbpara
            knpar = tblp(1+4*(j-1))
            nomjv = tblp(1+4*(j-1)+2)
            nomjvl = tblp(1+4*(j-1)+3)
            if (inpar .eq. knpar) then
                npara = npara+1
                nom_para(npara) = j
                call jeveuo(nomjv, 'L', val_para(npara))
                call jeveuo(nomjvl, 'L', log_para(npara))
                ilon = lxlgut(inpar)
                ilmp = max(ilon, ilmp)
                goto 10
            end if
        end do
        erreur = .true.
        valk = inpar
        call utmess('A', 'UTILITAI6_89', sk=valk)
10      continue
    end do
    if (erreur) then
        call utmess('F', 'PREPOST_60')
    end if
!
    nparaf = npara
    chaine = ' '
    chain2 = ' '
    ideb = 2
    do i = 1, npara
        ipar = nom_para(i)
        type = tblp(1+4*(ipar-1)+1)
        ilon = lxlgut(tblp(1+4*(ipar-1)))
        if (type(1:3) .eq. 'K80') then
            iaj = (80-ilon)/2
            ifin = ideb+80-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
            if (format .eq. 'AGRAF') ifin = ifin+1
        else if (type(1:1) .eq. 'I') then
            ilm = max(12, ilmp)
            iaj = (ilm-ilon)/2
            ifin = ideb+ilm-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
        else if (type(1:1) .eq. 'R') then
            ilm = max(ir, ilmp)
            iaj = (ilm-ilon)/2
            ifin = ideb+ilm-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
        else if (type(1:1) .eq. 'C') then
            ilm = 2*ir+1
            ilm = max(ilm, ilmp)
            iaj = (ilm-ilon)/2
            ifin = ideb+ilm-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj+8:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
        else if (type(1:2) .eq. 'K8') then
            ilm = max(8, ilmp)
            iaj = (ilm-ilon)/2
            ifin = ideb+ilm-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
            if (format .eq. 'AGRAF') ifin = ifin+1
        else if (type(1:3) .eq. 'K16') then
            ilm = max(16, ilmp)
            iaj = (ilm-ilon)/2
            ifin = ideb+ilm-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
            if (format .eq. 'AGRAF') ifin = ifin+1
        else if (type(1:3) .eq. 'K24') then
            iaj = (24-ilon)/2
            ifin = ideb+24-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
            if (format .eq. 'AGRAF') ifin = ifin+1
        else if (type(1:3) .eq. 'K32') then
            iaj = (32-ilon)/2
            ifin = ideb+32-1
            if (ifin .gt. 1999) then
                ifin = ideb
                nparaf = i-1
                goto 22
            end if
            chaine(ideb+iaj:ifin) = tblp(1+4*(ipar-1))
            chain2(ideb+iaj:ifin) = type
            if (format .eq. 'AGRAF') ifin = ifin+1
        end if
!
        ideb = ifin+2
    end do
22  continue
    if (nparaf .ne. npara) then
        call utmess('A', 'UTILITAI4_84')
    end if
    call codent(ifin, 'G', chfin)
    form1 = '(A'//chfin//')'
    write (ifr, form1) chaine(1:ifin)
!
    if (format .eq. 'ASTER') then
        write (ifr, form1) chain2(1:ifin)
    end if
!
    do i = 1, nblign
        chaine = ' '
        ideb = 2
        do j = 1, nparaf
            ipar = nom_para(j)
            type = tblp(1+4*(ipar-1)+1)
            jvale = val_para(j)
            jlogq = log_para(j)
            if (zi(jlogq+i-1) .eq. 1) then
                if (type(1:1) .eq. 'I') then
                    ilm = max(ilmp, 12)
                    ifin = ideb+ilm-1
                    write (chaine(ideb:ifin), '(I12)') zi(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:1) .eq. 'R') then
                    ilm = max(ilmp, ir)
                    ifin = ideb+ilm-1
                    write (chaine(ideb:ifin), formr) zr(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:1) .eq. 'C') then
                    ilm = 2*ir+1
                    ilm = max(ilm, ilmp)/2
                    ifin = ideb+ilm-1
                    write (chaine(ideb:ifin), formr) zc(jvale+i-1)
                    ifin = ideb+ilm-1
                    write (chaine(ideb:ifin), formr) zc(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:3) .eq. 'K80') then
                    if (format .eq. 'AGRAF') then
                        chaine(ideb:ideb) = bacs
                        ideb = ideb+1
                    end if
                    ifin = ideb+80-1
                    chaine(ideb:ifin) = zk80(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:3) .eq. 'K32') then
                    if (format .eq. 'AGRAF') then
                        chaine(ideb:ideb) = bacs
                        ideb = ideb+1
                    end if
                    ifin = ideb+32-1
                    chaine(ideb:ifin) = zk32(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:3) .eq. 'K24') then
                    if (format .eq. 'AGRAF') then
                        chaine(ideb:ideb) = bacs
                        ideb = ideb+1
                    end if
                    ifin = ideb+24-1
                    chaine(ideb:ifin) = zk24(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:3) .eq. 'K16') then
                    if (format .eq. 'AGRAF') then
                        chaine(ideb:ideb) = bacs
                        ideb = ideb+1
                    end if
                    ilm = max(ilmp, 16)
                    ifin = ideb+ilm-1
                    chaine(ideb:ifin) = zk16(jvale+i-1)
                    ideb = ifin+2
                else if (type(1:2) .eq. 'K8') then
                    if (format .eq. 'AGRAF') then
                        chaine(ideb:ideb) = bacs
                        ideb = ideb+1
                    end if
                    ilm = max(ilmp, 8)
                    ifin = ideb+ilm-1
                    chaine(ideb:ifin) = zk8(jvale+i-1)
                    ideb = ifin+2
                end if
            else
                if (type(1:3) .eq. 'K80') then
                    ifin = ideb+80-1
                    ideb = ideb+39
                else if (type(1:1) .eq. 'I') then
                    ilm = max(ilmp, 12)
                    ifin = ideb+ilm-1
                    ideb = ideb+5
                else if (type(1:1) .eq. 'R') then
                    ilm = max(ilmp, ir)
                    ifin = ideb+ilm-1
                    ideb = ideb+5
                else if (type(1:1) .eq. 'C') then
                    ifin = ideb+25-1
                    ideb = ideb+11
                else if (type(1:2) .eq. 'K8') then
                    ilm = max(ilmp, 8)
                    ifin = ideb+ilm-1
                    ideb = ideb+5
                else if (type(1:3) .eq. 'K16') then
                    ilm = max(ilmp, 16)
                    ifin = ideb+ilm-1
                    ideb = ideb+7
                else if (type(1:3) .eq. 'K24') then
                    ifin = ideb+24-1
                    ideb = ideb+11
                else if (type(1:3) .eq. 'K32') then
                    ifin = ideb+32-1
                    ideb = ideb+15
                end if
                if (format .eq. 'AGRAF') then
                    ideb = ideb-1
                    chaine(ideb:ideb+1) = bacs//'-'
                else
                    chaine(ideb:ideb) = '-'
                end if
                ideb = ifin+2
            end if
        end do
        call codent(ifin, 'G', chfin)
        form1 = '(A'//chfin//')'
        write (ifr, form1) chaine(1:ifin)
    end do
!
    AS_DEALLOCATE(vi=nom_para)
    AS_DEALLOCATE(vi=val_para)
    AS_DEALLOCATE(vi=log_para)
!
    call jedema()
!
!
end subroutine
