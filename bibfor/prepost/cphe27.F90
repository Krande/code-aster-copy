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

subroutine cphe27(maout, inc, jcnnpa, conloc, &
                  limane, jmacou, jmacsu, macou, &
                  macsu, ind, ind1)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"

!
    character(len=8), intent(in) :: maout
    integer(kind=8), intent(in) :: inc
    integer(kind=8), intent(in) :: jcnnpa
    character(len=24), intent(in) :: conloc
    character(len=24), intent(in) :: limane
    integer(kind=8), intent(in) :: jmacou
    integer(kind=8), intent(in) :: jmacsu
    integer(kind=8), intent(in) :: macou
    integer(kind=8), intent(in) :: macsu
    integer(kind=8), intent(out) :: ind
    integer(kind=8), intent(out) :: ind1
! -------------------------------------------------------------------------------------------------
!        CREATION DES NOUVEAUS NOUEDS ET NOUVELLE MAILLE CAS HEXA27/QUAD9
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
    integer(kind=8) :: patch, ino
    integer(kind=8) :: jlimane, jconloc
! -------------------------------------------------------------------------------------------------
    call jemarq()
!
    call jecroc(jexnum(maout//'.PATCH', inc+1))
    call jeecra(jexnum(maout//'.PATCH', inc+1), 'LONMAX', ival=2)
    call jeecra(jexnum(maout//'.PATCH', inc+1), 'LONUTI', ival=2)
    call jeveuo(jexnum(maout//'.PATCH', inc+1), 'E', patch)
! ----- Type de maille du patch ------------------------------------------------------------------
    zi(patch-1+1) = 27
! ----- DDL interne ------------------------------------------------------------------------------
    zi(patch-1+2) = zi(jmacou-1+9)
! ----- .CONOPA ----------------------------------------------------------------------------------
    zi(jcnnpa+zi(jmacou-1+9)-1) = inc
! --- NOUVEAUX ELEMENTS DE PEAU
    call jeecra(jexnum(conloc, ind), 'LONMAX', ival=9)
    call jeecra(jexnum(conloc, ind), 'LONUTI', ival=9)
    call jeveuo(jexnum(conloc, ind), 'E', jconloc)
    do ino = 1, 9
        zi(jconloc+ino-1) = zi(jmacou+ino-1)
    end do

! --- NOUVEAUX ELEMENTS DE CORPS
    call jeecra(jexnum(conloc, ind+1), 'LONMAX', ival=27)
    call jeecra(jexnum(conloc, ind+1), 'LONUTI', ival=27)
    call jeveuo(jexnum(conloc, ind+1), 'E', jconloc)
    do ino = 1, 27
        zi(jconloc+ino-1) = zi(jmacsu+ino-1)
    end do
! --- CONNECTIVITE ANCIENS NOUVEAUX ELEMENTS (Peau)

    call jeveuo(jexnum(limane, macou), 'E', jlimane)
    zi(jlimane+1-1) = ind
! --- INFO PATCH LIE
    zi(jlimane+2-1) = inc
! --- CONNECTIVITE ANCIENS NOUVEAUX ELEMENTS (Volume)

    call jeveuo(jexnum(limane, macsu), 'E', jlimane)
    zi(jlimane+1-1) = ind+1
! --- Nettoyage / mis à jour
    ind = ind+2
    ind1 = ind1+0
!
    call jedema()
end subroutine
