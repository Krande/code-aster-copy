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

subroutine cpmpq8_2(conloc, numa, indno, indma)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"

!
    integer(kind=8), intent(in) :: indma
    integer(kind=8), intent(in) :: indno
    integer(kind=8), intent(in) :: numa
    character(len=24), intent(in) :: conloc
!
!
! ----------------------------------------------------------------------
!         CREATION DES MAILLES DES NOUVELLES MAILLES DE PEAU
!         SUR LA FACE DE LA ZONE DE CONTACT ESCLAVE
!         CAS QUAD 8
! ----------------------------------------------------------------------
! IN        CONLOC  CONNECTIVITE LOCALE
! IN        NUMA    NUMERO DE LA MAILLE COURANTE
! IN        INDNO   INDICE DU PREMIER NOEUD AJOUTE
! IN        INDMA   INDICE DE LA PREMIERE MAILLE AJOUTEE
! ----------------------------------------------------------------------
    integer(kind=8) :: jconloc
! ----------------------------------------------------------------------
    call jemarq()
! ----------------------------------------------------------------------
    call jeecra(jexnum(conloc, indma), 'LONMAX', ival=3)
    call jeecra(jexnum(conloc, indma), 'LONUTI', ival=3)
    call jeveuo(jexnum(conloc, indma), 'E', jconloc)
    zi(jconloc+1-1) = zi(numa+1-1)
    zi(jconloc+2-1) = zi(numa+2-1)
    zi(jconloc+3-1) = indno
    call jeecra(jexnum(conloc, indma+1), 'LONMAX', ival=3)
    call jeecra(jexnum(conloc, indma+1), 'LONUTI', ival=3)
    call jeveuo(jexnum(conloc, indma+1), 'E', jconloc)
    zi(jconloc+1-1) = zi(numa+2-1)
    zi(jconloc+2-1) = zi(numa+3-1)
    zi(jconloc+3-1) = indno
    call jeecra(jexnum(conloc, indma+2), 'LONMAX', ival=3)
    call jeecra(jexnum(conloc, indma+2), 'LONUTI', ival=3)
    call jeveuo(jexnum(conloc, indma+2), 'E', jconloc)
    zi(jconloc+1-1) = zi(numa+3-1)
    zi(jconloc+2-1) = zi(numa+4-1)
    zi(jconloc+3-1) = indno
    call jeecra(jexnum(conloc, indma+3), 'LONMAX', ival=3)
    call jeecra(jexnum(conloc, indma+3), 'LONUTI', ival=3)
    call jeveuo(jexnum(conloc, indma+3), 'E', jconloc)
    zi(jconloc+1-1) = zi(numa+4-1)
    zi(jconloc+2-1) = zi(numa+1-1)
    zi(jconloc+3-1) = indno
! ----------------------------------------------------------------------
    call jedema()
end subroutine
