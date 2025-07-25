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

subroutine cpmcpt15_2(conloc, jmacsu, indno, indma, conneo)
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
    integer(kind=8) :: lino(15), jconloc
! -------------------------------------------------------------------------------------------------
    call jemarq()
! -------------------------------------------------------------------------------------------------
    if (conneo(1) .ne. 0 .and. conneo(2) .ne. 0 .and. conneo(5) .ne. 0 .and. conneo(4) .ne. 0) then
! -------------------------------------------------------------------------------------------------
        lino(1) = 2
        lino(2) = 3
        lino(3) = 1
        lino(4) = 5
        lino(5) = 6
        lino(6) = 4

        lino(7) = 8
        lino(8) = 9
        lino(9) = 7
        lino(10) = 11
        lino(11) = 12
        lino(12) = 10
        lino(13) = 14
        lino(14) = 15
        lino(15) = 13
        !write(*,*) '2'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(1) .ne. 0 .and. conneo(3) .ne. 0 .and. conneo(6) .ne. 0 .and. &
            conneo(4) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 1
        lino(2) = 2
        lino(3) = 3
        lino(4) = 4
        lino(5) = 5
        lino(6) = 6

        lino(7) = 7
        lino(8) = 8
        lino(9) = 9
        lino(10) = 10
        lino(11) = 11
        lino(12) = 12
        lino(13) = 13
        lino(14) = 14
        lino(15) = 15
        !write(*,*) '1'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(3) .ne. 0 .and. conneo(2) .ne. 0 .and. conneo(5) .ne. 0 .and. &
            conneo(6) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 3
        lino(2) = 1
        lino(3) = 2
        lino(4) = 6
        lino(5) = 4
        lino(6) = 5

        lino(7) = 9
        lino(8) = 7
        lino(9) = 8
        lino(10) = 12
        lino(11) = 10
        lino(12) = 11
        lino(13) = 15
        lino(14) = 13
        lino(15) = 14
        !write(*,*) '3'

    else
        ASSERT(.false.)
    end if
! -------------------------------------------------------------------------------------------------
    call jeecra(jexnum(conloc, indma), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(1)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+5-1) = indno

    zi(jconloc+6-1) = zi(jmacsu+lino(10)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(13)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(11)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(7)-1)
    zi(jconloc+10-1) = indno+conneo(lino(1))
    zi(jconloc+11-1) = indno+conneo(lino(4))
    zi(jconloc+12-1) = indno+4+2
    zi(jconloc+13-1) = indno+4+1
    !write(*,*) "ELEMENT",1
    !do ind=1, 13
    !    write(*,*)zi(jconloc+ind-1)
    !enddo

    call jeecra(jexnum(conloc, indma+1), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma+1), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma+1), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(3)-1)
    zi(jconloc+5-1) = indno

    zi(jconloc+6-1) = zi(jmacsu+lino(11)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(14)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(12)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(8)-1)
    zi(jconloc+10-1) = indno+4+1
    zi(jconloc+11-1) = indno+4+2
    zi(jconloc+12-1) = indno+conneo(lino(6))
    zi(jconloc+13-1) = indno+conneo(lino(3))
    !write(*,*) "ELEMENT",2
    !do ind=1, 13
    !    write(*,*)zi(jconloc+ind-1)
    !enddo
    call jeecra(jexnum(conloc, indma+2), 'LONMAX', ival=10)
    call jeecra(jexnum(conloc, indma+2), 'LONUTI', ival=10)
    call jeveuo(jexnum(conloc, indma+2), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+4-1) = indno

    zi(jconloc+5-1) = zi(jmacsu+lino(14)-1)
    zi(jconloc+6-1) = zi(jmacsu+lino(13)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(15)-1)
    zi(jconloc+8-1) = indno+conneo(lino(6))
    zi(jconloc+9-1) = indno+4+2
    zi(jconloc+10-1) = indno+conneo(lino(4))
    !write(*,*) "ELEMENT",3
    !do ind=1, 10
    !    write(*,*)zi(jconloc+ind-1)
    !enddo
    call jeecra(jexnum(conloc, indma+3), 'LONMAX', ival=10)
    call jeecra(jexnum(conloc, indma+3), 'LONUTI', ival=10)
    call jeveuo(jexnum(conloc, indma+3), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(1)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(3)-1)
    zi(jconloc+4-1) = indno

    zi(jconloc+5-1) = zi(jmacsu+lino(7)-1)
    zi(jconloc+6-1) = zi(jmacsu+lino(8)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(9)-1)
    zi(jconloc+8-1) = indno+conneo(lino(1))
    zi(jconloc+9-1) = indno+4+1
    zi(jconloc+10-1) = indno+conneo(lino(3))
    !write(*,*) "ELEMENT",4
    !do ind=1, 10
    !    write(*,*)zi(jconloc+ind-1)
    !enddo
! -------------------------------------------------------------------------------------------------
    call jedema()
end subroutine
