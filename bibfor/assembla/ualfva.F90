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
subroutine ualfva(mataz, basz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    character(len=*) :: mataz, basz
!     creation de l'objet mataz.VALM a partir de l'objet mataz.UALF
!     l'objet .UALF doit contenir la matrice initiale non factorisee :
!       - on cree l'objet.valm
!       - on detruit .ualf (a la fin de la routine)
!       - on cree le stockage morse dans le nume_ddl s'il n'existe pas.
!
!     cette routine ne devrait etre utilisee que rarement :
!        lorsque la matr_asse a ete cree sous la forme .ualf pour des
!        raisons historiques.
!
! ---------------------------------------------------------------------
! in  jxvar k19 mataz     : nom d'une s.d. matr_asse
! in        k1  basz      : base de creation pour .VALM
!                  si basz=' ' on prend la meme base que celle de .ualf
!     remarque : on detruit l'objet .UALF
!
!-----------------------------------------------------------------------
!     VARIABLES LOCALES
    character(len=1) :: base, tyrc
    character(len=14) :: nu
    character(len=19) :: stomor, matas
    integer(kind=8) :: neq, nbloc, nblocm, iret
    integer(kind=8) :: jsmhc
    integer(kind=8) :: itbloc, ieq, ibloc, jualf, jvale, kterm, nbterm, ilig
    integer(kind=8) :: ismdi, ismdi0, ibloav, iscdi, kblocm, nblocl
    integer(kind=8), pointer :: smde(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    character(len=24), pointer :: refa(:) => null()
!   ------------------------------------------------------------------
!
!
    call jemarq()
    matas = mataz
    base = basz
    if (base .eq. ' ') call jelira(matas//'.UALF', 'CLAS', cval=base)
!
!   -- .VALM ne doit pas exister :
    call jeexin(matas//'.VALM', iret)
    ASSERT(iret .eq. 0)
!
    call jeveuo(matas//'.REFA', 'L', vk24=refa)
    nu = refa(2) (1:14)
    stomor = nu//'.SMOS'
!
!   -- On ne sait traiter que les matrices generalisees :
    ASSERT(refa(10) .eq. 'GENE')
!
    call jeveuo(stomor//'.SMDE', 'L', vi=smde)
    neq = smde(1)
    nbloc = smde(3)
    call jeveuo(stomor//'.SMDI', 'L', vi=smdi)
    call jeveuo(stomor//'.SMHC', 'L', jsmhc)
    itbloc = smde(2)
!
    call jelira(matas//'.UALF', 'NMAXOC', nblocl)
    ASSERT(nblocl .eq. nbloc .or. nblocl .eq. 2*nbloc)
    nblocm = 1
    if (nblocl .eq. 2*nbloc) nblocm = 2
!
!   -- reel ou complexe ?
    call jelira(matas//'.UALF', 'TYPE', cval=tyrc)
    ASSERT(tyrc .eq. 'R' .or. tyrc .eq. 'C')
!
!
!     1. Allocation de .VALM :
!     ----------------------------------------
    call jecrec(matas//'.VALM', base//' V '//tyrc, 'NU', 'DISPERSE', 'CONSTANT', &
                nblocm)
    call jeecra(matas//'.VALM', 'LONMAX', itbloc)
    do kblocm = 1, nblocm
        call jecroc(jexnum(matas//'.VALM', kblocm))
    end do
!
!
!     2. Remplissage de .VALM :
!     ----------------------------------------
    do kblocm = 1, nblocm
        call jeveuo(jexnum(matas//'.VALM', kblocm), 'E', jvale)
        ibloav = 0+nbloc*(kblocm-1)
        ismdi0 = 0
        do ieq = 1, neq
            iscdi = smdi(ieq)
            ibloc = 1+nbloc*(kblocm-1)
!
!          -- on ramene le bloc en memoire si necessaire:
            if (ibloc .ne. ibloav) then
                call jeveuo(jexnum(matas//'.UALF', ibloc), 'L', jualf)
                if (ibloav .ne. 0) then
                    call jelibe(jexnum(matas//'.UALF', ibloav))
                end if
                ibloav = ibloc
            end if
!
            ismdi = smdi(ieq)
            nbterm = ismdi-ismdi0
!
            do kterm = 1, nbterm
                ilig = zi4(jsmhc-1+ismdi0+kterm)
                if (tyrc .eq. 'R') then
                    zr(jvale-1+ismdi0+kterm) = zr(jualf-1+iscdi+ilig- &
                                                  ieq)
                else
                    zc(jvale-1+ismdi0+kterm) = zc(jualf-1+iscdi+ilig- &
                                                  ieq)
                end if
            end do
            ASSERT(ilig .eq. ieq)
!
            ismdi0 = ismdi
        end do
    end do
!
!
!
    call jedetr(matas//'.UALF')
!
    call jedema()
!     CALL CHEKSD('sd_matr_asse',MATAS,IRET)
end subroutine
