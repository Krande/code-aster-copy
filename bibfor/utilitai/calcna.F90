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
subroutine calcna(nomfin, nomfon, nbvalp, valep, noparp, &
                  nbvalf, valef, noparf)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbvalp, nbvalf
    real(kind=8) :: valep(*), valef(*)
    character(len=19) :: nomfin, nomfon
    character(len=24) :: noparp, noparf
!
!     CREATION DU SD FONCTION A PARTIR D'UNE FORMULE (NAPPE )
!     ------------------------------------------------------------------
    integer(kind=8) :: lont, i, ival, lval, lfon, lprol, lpara, ier
    real(kind=8) :: vale(2)
    character(len=16) :: nopara(2)
!     ------------------------------------------------------------------
!
    call jemarq()
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NOMFON.VALE ---
!
    lont = 2*nbvalf*nbvalp
    nopara(1) = noparf
    nopara(2) = noparp
!
    call jecrec(nomfon//'.VALE', ' G V R', 'NU', 'CONTIG', 'VARIABLE', &
                nbvalp)
    call jeecra(nomfon//'.VALE', 'LONT', lont)
    do i = 1, nbvalp
        call jecroc(jexnum(nomfon//'.VALE', i))
        call jeecra(jexnum(nomfon//'.VALE', i), 'LONMAX', 2*nbvalf)
        call jeecra(jexnum(nomfon//'.VALE', i), 'LONUTI', 2*nbvalf)
        call jeveuo(jexnum(nomfon//'.VALE', i), 'E', lval)
        lfon = lval+nbvalf
        vale(2) = valep(i)
        do ival = 0, nbvalf-1
            zr(lval+ival) = valef(ival+1)
            vale(1) = zr(lval+ival)
            call fointe('F', nomfin, 2, nopara, vale, &
                        zr(lfon+ival), ier)
        end do
    end do
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NOMFON.PROL ---
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 7+2*nbvalp, lprol)
!
    zk24(lprol) = 'NAPPE           '
    zk24(lprol+1) = 'LIN LIN         '
    zk24(lprol+2) = noparp
    zk24(lprol+3) = 'TOUTRESU        '
    zk24(lprol+4) = 'EE              '
    zk24(lprol+5) = nomfon
    zk24(lprol+6) = noparf
    do ival = 1, nbvalp
        zk24(lprol+6+(2*ival-1)) = 'LIN LIN         '
        zk24(lprol+6+(2*ival)) = 'EE              '
    end do
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NOMFON.PARA ---
!
    call wkvect(nomfon//'.PARA', 'G V R', nbvalp, lpara)
    do ival = 1, nbvalp
        zr(lpara+ival-1) = valep(ival)
    end do
!
    call jedema()
end subroutine
