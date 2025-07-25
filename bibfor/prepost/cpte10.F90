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

subroutine cpte10(main, maout, inc, jcoor, jcnnpa, conloc, &
                  limane, nomnoe, nbno, jmacou, jmacsu, macou, &
                  macsu, ind, ind1)
!
    implicit none
#include "jeveux.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/cpnoq.h"
#include "asterfort/cnpc.h"
#include "asterfort/cpnov.h"
!
    character(len=8), intent(in) :: main
    character(len=8), intent(in) :: maout
    integer(kind=8), intent(in) :: inc
    integer(kind=8), intent(in) :: jcoor
    integer(kind=8), intent(in) :: jcnnpa
    character(len=24), intent(in) :: conloc
    character(len=24), intent(in) :: limane
    character(len=24), intent(in) :: nomnoe
    integer(kind=8), intent(in) :: nbno
    integer(kind=8), intent(in) :: jmacou
    integer(kind=8), intent(in) :: jmacsu
    integer(kind=8), intent(in) :: macou
    integer(kind=8), intent(in) :: macsu
    integer(kind=8), intent(out) :: ind
    integer(kind=8), intent(out) :: ind1
! -------------------------------------------------------------------------------------------------
!        CREATION DES NOUVEAUS NOUEDS ET NOUVELLE MAILLE CAS TETRA 10
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
    integer(kind=8) :: patch
    integer(kind=8) :: jlimane, jconloc
    integer(kind=8) :: jconneo
    character(len=24) :: conneo
! -------------------------------------------------------------------------------------------------
    call jemarq()
    call jecroc(jexnum(maout//'.PATCH', inc+1))
    call jeecra(jexnum(maout//'.PATCH', inc+1), 'LONMAX', ival=2)
    call jeecra(jexnum(maout//'.PATCH', inc+1), 'LONUTI', ival=2)
    call jeveuo(jexnum(maout//'.PATCH', inc+1), 'E', patch)
! --- TYPE DE MAILLE PATCH
    zi(patch-1+1) = 19
! --- DDL INTERNE
    zi(patch-1+2) = nbno+ind1+3
    zi(jcnnpa+nbno+ind1+3-1) = inc
! --- CREATION DES NOEUDS DDL INTERNE
    call cpnoq(main, macou, zr(jcoor), nbno+ind1, nomnoe)
! --- NOUVEAUX ELEMENTS DE PEAU
    call jeecra(jexnum(conloc, ind), 'LONMAX', ival=6)
    call jeecra(jexnum(conloc, ind), 'LONUTI', ival=6)
    call jeveuo(jexnum(conloc, ind), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacou+1-1)
    zi(jconloc+2-1) = zi(jmacou+2-1)
    zi(jconloc+3-1) = nbno+ind1+3
    zi(jconloc+4-1) = zi(jmacou+4-1)
    zi(jconloc+5-1) = nbno+ind1+1
    zi(jconloc+6-1) = nbno+ind1
    call jeecra(jexnum(conloc, ind+1), 'LONMAX', ival=6)
    call jeecra(jexnum(conloc, ind+1), 'LONUTI', ival=6)
    call jeveuo(jexnum(conloc, ind+1), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacou+2-1)
    zi(jconloc+2-1) = zi(jmacou+3-1)
    zi(jconloc+3-1) = nbno+ind1+3
    zi(jconloc+4-1) = zi(jmacou+5-1)
    zi(jconloc+5-1) = nbno+ind1+2
    zi(jconloc+6-1) = nbno+ind1+1
    call jeecra(jexnum(conloc, ind+2), 'LONMAX', ival=6)
    call jeecra(jexnum(conloc, ind+2), 'LONUTI', ival=6)
    call jeveuo(jexnum(conloc, ind+2), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacou+3-1)
    zi(jconloc+2-1) = zi(jmacou+1-1)
    zi(jconloc+3-1) = nbno+ind1+3
    zi(jconloc+4-1) = zi(jmacou+6-1)
    zi(jconloc+5-1) = nbno+ind1
    zi(jconloc+6-1) = nbno+ind1+2
! --- NOUVEAUX ELEMENTS DE CORPS
    conneo = '&&CPTE10.CNORD'
    call cnpc(main, macou, macsu, conneo)
    call jeveuo(conneo, 'L', jconneo)
    call cpnov(main, macsu, zr(jcoor), nbno+ind1+4, nomnoe, zi(jconneo))
    if (zi(jconneo+1-1) .eq. 0) then
        call jeecra(jexnum(conloc, ind+3), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+3), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+3), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+1-1)
        zi(jconloc+2-1) = zi(jmacsu+2-1)
        zi(jconloc+3-1) = zi(jmacsu+3-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+5-1)
        zi(jconloc+6-1) = zi(jmacsu+6-1)
        zi(jconloc+7-1) = zi(jmacsu+7-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+2-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+3-1)-1
        call jeecra(jexnum(conloc, ind+4), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+4), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+4), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+1-1)
        zi(jconloc+2-1) = zi(jmacsu+3-1)
        zi(jconloc+3-1) = zi(jmacsu+4-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+7-1)
        zi(jconloc+6-1) = zi(jmacsu+10-1)
        zi(jconloc+7-1) = zi(jmacsu+8-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+3-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+4-1)-1
        call jeecra(jexnum(conloc, ind+5), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+5), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+5), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+1-1)
        zi(jconloc+2-1) = zi(jmacsu+4-1)
        zi(jconloc+3-1) = zi(jmacsu+2-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+8-1)
        zi(jconloc+6-1) = zi(jmacsu+9-1)
        zi(jconloc+7-1) = zi(jmacsu+5-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+4-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+2-1)-1
    elseif (zi(jconneo+2-1) .eq. 0) then
        call jeecra(jexnum(conloc, ind+3), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+3), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+3), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+2-1)
        zi(jconloc+2-1) = zi(jmacsu+3-1)
        zi(jconloc+3-1) = zi(jmacsu+1-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+6-1)
        zi(jconloc+6-1) = zi(jmacsu+7-1)
        zi(jconloc+7-1) = zi(jmacsu+5-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+3-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+1-1)-1
        call jeecra(jexnum(conloc, ind+4), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+4), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+4), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+2-1)
        zi(jconloc+2-1) = zi(jmacsu+4-1)
        zi(jconloc+3-1) = zi(jmacsu+3-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+9-1)
        zi(jconloc+6-1) = zi(jmacsu+10-1)
        zi(jconloc+7-1) = zi(jmacsu+6-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+4-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+3-1)-1
        call jeecra(jexnum(conloc, ind+5), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+5), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+5), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+2-1)
        zi(jconloc+2-1) = zi(jmacsu+1-1)
        zi(jconloc+3-1) = zi(jmacsu+4-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+5-1)
        zi(jconloc+6-1) = zi(jmacsu+8-1)
        zi(jconloc+7-1) = zi(jmacsu+9-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+1-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+4-1)-1
    elseif (zi(jconneo+3-1) .eq. 0) then
        call jeecra(jexnum(conloc, ind+3), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+3), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+3), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+3-1)
        zi(jconloc+2-1) = zi(jmacsu+1-1)
        zi(jconloc+3-1) = zi(jmacsu+2-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+7-1)
        zi(jconloc+6-1) = zi(jmacsu+5-1)
        zi(jconloc+7-1) = zi(jmacsu+6-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+1-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+2-1)-1
        call jeecra(jexnum(conloc, ind+4), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+4), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+4), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+3-1)
        zi(jconloc+2-1) = zi(jmacsu+4-1)
        zi(jconloc+3-1) = zi(jmacsu+1-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+10-1)
        zi(jconloc+6-1) = zi(jmacsu+8-1)
        zi(jconloc+7-1) = zi(jmacsu+7-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+4-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+1-1)-1
        call jeecra(jexnum(conloc, ind+5), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+5), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+5), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+3-1)
        zi(jconloc+2-1) = zi(jmacsu+2-1)
        zi(jconloc+3-1) = zi(jmacsu+4-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+6-1)
        zi(jconloc+6-1) = zi(jmacsu+9-1)
        zi(jconloc+7-1) = zi(jmacsu+10-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+2-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+4-1)-1
    elseif (zi(jconneo+4-1) .eq. 0) then
        call jeecra(jexnum(conloc, ind+3), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+3), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+3), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+4-1)
        zi(jconloc+2-1) = zi(jmacsu+1-1)
        zi(jconloc+3-1) = zi(jmacsu+3-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+8-1)
        zi(jconloc+6-1) = zi(jmacsu+7-1)
        zi(jconloc+7-1) = zi(jmacsu+10-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+1-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+3-1)-1
        call jeecra(jexnum(conloc, ind+4), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+4), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+4), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+4-1)
        zi(jconloc+2-1) = zi(jmacsu+2-1)
        zi(jconloc+3-1) = zi(jmacsu+1-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+9-1)
        zi(jconloc+6-1) = zi(jmacsu+5-1)
        zi(jconloc+7-1) = zi(jmacsu+8-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+2-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+1-1)-1
        call jeecra(jexnum(conloc, ind+5), 'LONMAX', ival=10)
        call jeecra(jexnum(conloc, ind+5), 'LONUTI', ival=10)
        call jeveuo(jexnum(conloc, ind+5), 'E', jconloc)
        zi(jconloc+1-1) = zi(jmacsu+4-1)
        zi(jconloc+2-1) = zi(jmacsu+3-1)
        zi(jconloc+3-1) = zi(jmacsu+2-1)
        zi(jconloc+4-1) = nbno+ind1+3
        zi(jconloc+5-1) = zi(jmacsu+10-1)
        zi(jconloc+6-1) = zi(jmacsu+6-1)
        zi(jconloc+7-1) = zi(jmacsu+9-1)
        zi(jconloc+8-1) = nbno+ind1+4
        zi(jconloc+9-1) = nbno+ind1+zi(jconneo+3-1)-1
        zi(jconloc+10-1) = nbno+ind1+zi(jconneo+2-1)-1
    end if
! --- CONNECTIVITE ANCIENS NOUVEAUX ELEMENTS

    call jeveuo(jexnum(limane, macou), 'E', jlimane)
    zi(jlimane+1-1) = ind
    zi(jlimane+2-1) = ind+1
    zi(jlimane+3-1) = ind+2
! --- INFO PATCH LIE
    zi(jlimane+4-1) = inc
!

    call jeveuo(jexnum(limane, macsu), 'E', jlimane)
    zi(jlimane+1-1) = ind+3
    zi(jlimane+2-1) = ind+4
    zi(jlimane+3-1) = ind+5
    ind = ind+6
    ind1 = ind1+5
    call jedetr(conneo)
!
    call jedema()
end subroutine
