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

subroutine cphe20(main, maout, inc, jcoor, jcnnpa, conloc, &
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
#include "asterfort/cpnph20.h"
#include "asterfort/cpnch20.h"
#include "asterfort/cpmph20.h"
#include "asterfort/cpmch20.h"
#include "asterfort/cnpc.h"
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
    integer(kind=8) :: jlimane
    integer(kind=8) :: jconneo
    character(len=24) :: conneo
! -------------------------------------------------------------------------------------------------
    call jemarq()
!
    call jecroc(jexnum(maout//'.PATCH', inc+1))
    call jeecra(jexnum(maout//'.PATCH', inc+1), 'LONMAX', ival=5)
    call jeecra(jexnum(maout//'.PATCH', inc+1), 'LONUTI', ival=5)
    call jeveuo(jexnum(maout//'.PATCH', inc+1), 'E', patch)
! --- TYPE DE MAILLE PATCH
    zi(patch-1+1) = 26
! --- DDL INTERNE
    zi(patch-1+2) = nbno+ind1
    zi(jcnnpa+nbno+ind1-1) = inc
! --- DDLs SUPPLEMENTAIRES
    zi(patch-1+3) = nbno+ind1+1
    zi(jcnnpa+nbno+ind1+1-1) = inc
    zi(patch-1+4) = nbno+ind1+2
    zi(jcnnpa+nbno+ind1+2-1) = inc
    zi(patch-1+5) = nbno+ind1+3
    zi(jcnnpa+nbno+ind1+3-1) = inc
! --- CREATION DES NOEUDS DDL INTERNE
    call cpnph20(main, macou, zr(jcoor), nbno+ind1, nomnoe)
! --- NOUVEAUX ELEMENTS DE PEAU
    call cpmph20(conloc, jmacou, nbno+ind1, ind)
! --- CREATION DES NOEUDS DDL DANS LE VOLUME
    conneo = '&&CPHE20.CNORD'
    call cnpc(main, macou, macsu, conneo)
    call jeveuo(conneo, 'L', jconneo)
    call cpnch20(main, macsu, zr(jcoor), nbno+ind1+12, nomnoe, zi(jconneo))
! --- NOUVEAUX ELEMENTS DE CORPS
    call cpmch20(conloc, jmacsu, nbno+ind1, ind+5, zi(jconneo))
! --- CONNECTIVITE ANCIENS NOUVEAUX ELEMENTS (Peau)

    call jeveuo(jexnum(limane, macou), 'E', jlimane)
    zi(jlimane+1-1) = ind
    zi(jlimane+2-1) = ind+1
    zi(jlimane+3-1) = ind+2
    zi(jlimane+4-1) = ind+3
    zi(jlimane+5-1) = ind+4
! --- INFO PATCH LIE
    zi(jlimane+6-1) = inc
! --- CONNECTIVITE ANCIENS NOUVEAUX ELEMENTS (Volume)

    call jeveuo(jexnum(limane, macsu), 'E', jlimane)
    zi(jlimane+1-1) = ind+5
    zi(jlimane+2-1) = ind+6
    zi(jlimane+3-1) = ind+7
    zi(jlimane+4-1) = ind+8
    zi(jlimane+5-1) = ind+9
    zi(jlimane+6-1) = ind+10
! --- Nettoyage / mis à jour
    ind = ind+11
    ind1 = ind1+28
    call jedetr(conneo)
!
    call jedema()
end subroutine
