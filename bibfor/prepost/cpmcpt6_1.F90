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

subroutine cpmcpt6_1(conloc, jmacsu, indno, indma, conneo)
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
#include "asterfort/assert.h"

!
    integer(kind=8), intent(in) :: indma
    integer(kind=8), intent(in) :: indno
    integer(kind=8), intent(in) :: jmacsu
    integer(kind=8), intent(in) :: conneo(*)
    character(len=24), intent(in) :: conloc
!
! -------------------------------------------------------------------------------------------------
!         CREATION DES MAILLES DES NOUVELLES MAILLES DE PEAU
!         SUR LA FACE DE LA ZONE DE CONTACT ESCLAVE
!         CAS QUAD 8
! -------------------------------------------------------------------------------------------------
! IN        CONLOC  CONNECTIVITE LOCALE
! IN        NUMA    NUMERO DE LA MAILLE COURANTE
! IN        INDNO   INDICE DU PREMIER NOEUD AJOUTE
! IN        INDMA   INDICE DE LA PREMIERE MAILLE AJOUTEE
! -------------------------------------------------------------------------------------------------
    integer(kind=8) :: lino(8), jconloc
! -------------------------------------------------------------------------------------------------
    call jemarq()
! -------------------------------------------------------------------------------------------------
    if (conneo(1) .ne. 0 .and. conneo(2) .ne. 0 .and. conneo(3) .ne. 0) then
! -------------------------------------------------------------------------------------------------
        lino(1) = 1
        lino(2) = 2
        lino(3) = 3
        lino(4) = 4
        lino(5) = 5
        lino(6) = 6
        !write(*,*) '1'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(4) .ne. 0 .and. conneo(5) .ne. 0 .and. conneo(6) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 6
        lino(2) = 5
        lino(3) = 4
        lino(4) = 3
        lino(5) = 2
        lino(6) = 1
        !write(*,*) '2'

    else
        ASSERT(.false.)
    end if
! -------------------------------------------------------------------------------------------------
    call jeecra(jexnum(conloc, indma), 'LONMAX', ival=5)
    call jeecra(jexnum(conloc, indma), 'LONUTI', ival=5)
    call jeveuo(jexnum(conloc, indma), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(1)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+5-1) = indno
    call jeecra(jexnum(conloc, indma+1), 'LONMAX', ival=5)
    call jeecra(jexnum(conloc, indma+1), 'LONUTI', ival=5)
    call jeveuo(jexnum(conloc, indma+1), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(1)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(3)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+5-1) = indno
    call jeecra(jexnum(conloc, indma+2), 'LONMAX', ival=5)
    call jeecra(jexnum(conloc, indma+2), 'LONUTI', ival=5)
    call jeveuo(jexnum(conloc, indma+2), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(3)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+5-1) = indno
    call jeecra(jexnum(conloc, indma+3), 'LONMAX', ival=4)
    call jeecra(jexnum(conloc, indma+3), 'LONUTI', ival=4)
    call jeveuo(jexnum(conloc, indma+3), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+4-1) = indno
! -------------------------------------------------------------------------------------------------
    call jedema()
end subroutine
