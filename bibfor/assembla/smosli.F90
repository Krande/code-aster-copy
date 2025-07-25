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

subroutine smosli(stomoz, stolcz, basz, rtbloc)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/ldlt_renum.h"
    character(len=*) :: stomoz, stolcz, basz
    real(kind=8) :: rtbloc
!     calcul d'un stoc_lciel a partir d'un stoc_morse (pour contenir
!     la meme matrice)
!     ------------------------------------------------------------------
! in  jxin  k19 stomoz     : nom d'une s.d. stoc_morse
! in  jxout k19 stolcz     : nom d'une s.d. stoc_lciel
! in        k1  basz       : base de creation pour stolcz
! in        r   rtbloc     : taille des blocs de stolci
!-----------------------------------------------------------------------

    character(len=1) :: base
    character(len=19) :: stomor, stolci
    integer(kind=8) ::  neq, nbloc, hcc, nnz, kterm, jcolm, jcollc, iligm, iliglc, temp
    integer(kind=8) :: hc, hcmax, itbloc, ieqlc, ibloc, tcumu
    integer(kind=8), pointer :: schc(:) => null()
    integer(kind=8), pointer :: scdi(:) => null()
    integer(kind=8), pointer :: scib(:) => null()
    integer(kind=8), pointer :: scbl(:) => null()
    integer(kind=8), pointer :: scde(:) => null()

    integer(kind=8), pointer :: smde(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=4), pointer :: smhc(:) => null()

    integer(kind=8), pointer :: m2lc(:) => null()
    integer(kind=8), pointer :: lc2m(:) => null()
!   ------------------------------------------------------------------
!
!
!
    call jemarq()
    stomor = stomoz
    stolci = stolcz
    ASSERT(stomor(1:14) .eq. stolci(1:14))
    base = basz
!
!   -- on detruit stolci s'il existe deja :
    call detrsd('STOC_LCIEL', stolci)
!
!
!   -- objet .SCDE : c'est facile mais pas complet:
!   ------------------------------------------------
    call wkvect(stolci//'.SCDE', base//' V I', 6, vi=scde)
    call jeveuo(stomor//'.SMDE', 'L', vi=smde)
    neq = smde(1)
    nnz = smde(2)
    scde(1) = neq

!   -- calcul de itbloc :
    itbloc = nint(rtbloc*1024)

!   -- calcul des vecteurs d'indirection :
!       .M2LC : ieqm  -> ieqlc
!       .LC2M : ieqlc -> ieqm
!   -----------------------------------------------
    call ldlt_renum(stomor(1:14))
    call jeveuo(stolci//'.M2LC', 'L', vi=m2lc)
    call jeveuo(stolci//'.LC2M', 'L', vi=lc2m)
    ASSERT(m2lc(1) .ge. 1 .and. m2lc(1) .le. neq)

!   -- allocation de  .SCHC .SCDI et .SCIB :
!   ------------------------------------------
    call wkvect(stolci//'.SCHC', base//' V I', neq, vi=schc)
    call wkvect(stolci//'.SCDI', base//' V I', neq, vi=scdi)
    call wkvect(stolci//'.SCIB', base//' V I', neq, vi=scib)
!
!
!   1. remplissage de .SCHC .SCDI et .SCIB
!   -------------------------------------------------------------
    call jeveuo(stomor//'.SMDI', 'L', vi=smdi)
    call jeveuo(stomor//'.SMHC', 'L', vi4=smhc)

!
!   -- 1.1  Calcul des hauteurs de colonnes  .SCHC :
!   -------------------------------------------------
    jcolm = 1
    do kterm = 1, nnz
        if (smdi(jcolm) .lt. kterm) jcolm = jcolm+1
        ASSERT(jcolm .ge. 0 .and. jcolm .le. neq)
        iligm = smhc(kterm)
        ASSERT(iligm .ge. 0 .and. iligm .le. neq)
        iliglc = m2lc(iligm)
        jcollc = m2lc(jcolm)
        if (iliglc .gt. jcollc) then
            temp = iliglc
            iliglc = jcollc
            jcollc = temp
        end if
        hc = jcollc+1-iliglc
        schc(jcollc) = max(schc(jcollc), hc)
    end do
    ASSERT(schc(1) .eq. 1)

!
!   -- 1.2  Calcul de .SCIB, .SCDI, hcmax et hcc :
!   ---------------------------------------------------

    hcmax = 0
    tcumu = 0
    hcc = 0
    nbloc = 1
    do ieqlc = 1, neq
        hc = schc(ieqlc)
!
!       -- peut-on encore stocker cette colonne dans le bloc courant ?
        ASSERT(hc .le. itbloc)
        if (tcumu+hc .gt. itbloc) then
            nbloc = nbloc+1
            tcumu = 0
        end if
!
        scib(ieqlc) = nbloc
        scdi(ieqlc) = tcumu+hc
        hcmax = max(hcmax, hc)
        tcumu = tcumu+hc
        hcc = hcc+hc
    end do
    scde(3) = nbloc
    scde(4) = hcmax
!
!
!   2. allocation et remplissage de .SCBL :
!   ----------------------------------------
    call wkvect(stolci//'.SCBL', base//' V I', nbloc+1, vi=scbl)
    scbl(1) = 0
    ibloc = 1
    do ieqlc = 1, neq
        if (scib(ieqlc) .gt. ibloc) then
            ibloc = ibloc+1
            scbl(ibloc) = ieqlc-1
        end if
    end do
    ASSERT(ibloc .eq. nbloc)
    scbl(nbloc+1) = neq
!
!
!   3. si toute la matrice tient dans un seul bloc, on
!      adapte la taille de ce bloc
!   ------------------------------------------------------
    if (itbloc .gt. hcc) then
        ASSERT(nbloc .eq. 1)
        itbloc = hcc
    end if
    scde(2) = itbloc
!
!
    call jedema()
end subroutine
