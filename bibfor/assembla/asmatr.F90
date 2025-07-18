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
subroutine asmatr(nbmat, tlimat, licoef, nu, &
                  listLoadZ, cumul, base, itysca, mataz)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/ascima.h"
#include "asterfort/assert.h"
#include "asterfort/assmam.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedbg2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/masyns.h"
#include "asterfort/typmat.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbmat, itysca
    character(len=*) :: base, mataz, licoef, nu
    character(len=19) :: tlimat(nbmat)
    character(len=*) :: listLoadZ
    character(len=4) :: cumul
!
! --------------------------------------------------------------------------------------------------
!
! in  i   nbmat  : nombre de matr_elem de la liste tlimat
! in  k19 tlimat : liste des matr_elem
! in  k24 licoef : nom du vecteur contenant les coef. mult.
!                  des matr_elem
!                  si licoef=' ' on prend 1.d0 comme coef.
! in  k14 nu     : nom du nume_ddl
! in  k19 listLoad : pour les charges cinematiques :
!                  / sd_infcha (k19)
!                  / nom d'un objet jeveux (k24) contenant
!                    les noms des charges cinematiques (k24)
! in  k4 cumul  : 'ZERO' ou 'CUMU'
!                 'ZERO':si un objet de nom matas et de type
!                        matr_asse existe on ecrase son contenu.
!                 'CUMU':si un objet de nom matas et de type
!                        matr_asse existe on cumule dans .valm
! in  k1  base   : base sur laquelle on cree l'objet mataz
! in  i   itysca  : type des matrices elementaires a assembler
!                          1 --> reelles
!                          2 --> complexes
! in/out k19 mataz : l'objet mataz de type matr_asse est cree et rempli
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1) :: matsym
    character(len=24) :: licoe2
    integer(kind=8) :: k
    character(len=19) :: tlima2(150), matas, infc19
    integer(kind=8) :: ilicoe, i, iret, ibid, idbgav
    integer(kind=8) :: jrefa
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call jedbg2(idbgav, 0)

    matas = mataz
    licoe2 = licoef
    infc19 = listLoadZ

    ASSERT(cumul .eq. 'ZERO' .or. cumul .eq. 'CUMU')
    if (cumul .eq. 'ZERO') call detrsd('MATR_ASSE', matas)

    ASSERT(nbmat .le. 150)
    do k = 1, nbmat
        tlima2(k) = tlimat(k)
    end do

!   -- traitement de la liste des coef. multiplicateurs :
!   ---------------------------------------------------------------
    if (licoe2 .eq. ' ') then
        call wkvect('&&ASMATR.LICOEF', 'V V R', nbmat, ilicoe)
        do i = 1, nbmat
            zr(ilicoe+i-1) = 1.d0
        end do
    else
        call jeveuo(licoe2, 'L', ilicoe)
    end if

!   -- preparation de la liste de matr_elem pour qu'ils soient
!      du meme type (symetrique ou non) que la matr_asse :
!   ---------------------------------------------------------------
    matsym = typmat(nbmat, tlima2)

!   -- si matrice existe deja et qu'elle doit etre non-symetrique,
!      on la de-symetrise :
!   ---------------------------------------------------------------
    if (cumul .eq. 'CUMU') then
        call jeexin(matas//'.REFA', iret)
        ASSERT(iret .gt. 0)
        call jeveuo(matas//'.REFA', 'L', jrefa)
        if (matsym .eq. 'N' .and. zk24(jrefa-1+9) .eq. 'MS') call masyns(matas)
    end if

!   -- assemblage proprement dit :
!   -------------------------------
    call assmam(base, matas, nbmat, tlima2, zr(ilicoe), &
                nu, cumul, itysca)

!   -- traitement des charges cinematiques :
!   ----------------------------------------
    call jeveuo(matas//'.REFA', 'L', jrefa)
    ASSERT(zk24(jrefa-1+3) .ne. 'ELIMF')
    call ascima(infc19, nu, matas, cumul)

!   -- menage :
!   -----------
    call jedetr('&&ASMATR.LICOEF')
    call jedetr(matas//'.&INT')
    call jedetr(matas//'.&IN2')
    call jedetr(matas//'.LILI')

    call jedbg2(ibid, idbgav)
    call jedema()
end subroutine
