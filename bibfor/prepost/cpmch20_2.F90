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

subroutine cpmch20_2(conloc, jmacsu, indno, indma, conneo)
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
    integer(kind=8) :: lino(20), jconloc
! -------------------------------------------------------------------------------------------------
    call jemarq()
! -------------------------------------------------------------------------------------------------
    if (conneo(1) .ne. 0 .and. conneo(2) .ne. 0 .and. conneo(3) .ne. 0 .and. conneo(4) .ne. 0) then
! -------------------------------------------------------------------------------------------------
        lino(1) = 4
        lino(2) = 3
        lino(3) = 2
        lino(4) = 1
        lino(5) = 8
        lino(6) = 7
        lino(7) = 6
        lino(8) = 5
        lino(9) = 11
        lino(10) = 10
        lino(11) = 9
        lino(12) = 12
        lino(13) = 16
        lino(14) = 15
        lino(15) = 14
        lino(16) = 13
        lino(17) = 19
        lino(18) = 18
        lino(19) = 17
        lino(20) = 20
        !write(*,*) '1'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(1) .ne. 0 .and. conneo(2) .ne. 0 .and. &
            conneo(6) .ne. 0 .and. conneo(5) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 1
        lino(2) = 2
        lino(3) = 6
        lino(4) = 5
        lino(5) = 4
        lino(6) = 3
        lino(7) = 7
        lino(8) = 8
        lino(9) = 9
        lino(10) = 14
        lino(11) = 17
        lino(12) = 13
        lino(13) = 12
        lino(14) = 10
        lino(15) = 18
        lino(16) = 20
        lino(17) = 11
        lino(18) = 15
        lino(19) = 19
        lino(20) = 16
        !write(*,*) '2'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(2) .ne. 0 .and. conneo(3) .ne. 0 .and. &
            conneo(7) .ne. 0 .and. conneo(6) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 2
        lino(2) = 3
        lino(3) = 7
        lino(4) = 6
        lino(5) = 1
        lino(6) = 4
        lino(7) = 8
        lino(8) = 5
        lino(9) = 10
        lino(10) = 15
        lino(11) = 18
        lino(12) = 14
        lino(13) = 9
        lino(14) = 11
        lino(15) = 19
        lino(16) = 17
        lino(17) = 12
        lino(18) = 16
        lino(19) = 20
        lino(20) = 13
        !write(*,*) '3'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(3) .ne. 0 .and. conneo(4) .ne. 0 .and. &
            conneo(8) .ne. 0 .and. conneo(7) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 3
        lino(2) = 4
        lino(3) = 8
        lino(4) = 7
        lino(5) = 2
        lino(6) = 1
        lino(7) = 5
        lino(8) = 6
        lino(9) = 11
        lino(10) = 16
        lino(11) = 19
        lino(12) = 15
        lino(13) = 10
        lino(14) = 12
        lino(15) = 20
        lino(16) = 18
        lino(17) = 9
        lino(18) = 13
        lino(19) = 17
        lino(20) = 14
        !write(*,*) '4'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(1) .ne. 0 .and. conneo(5) .ne. 0 .and. &
            conneo(8) .ne. 0 .and. conneo(4) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 1
        lino(2) = 5
        lino(3) = 8
        lino(4) = 4
        lino(5) = 2
        lino(6) = 6
        lino(7) = 7
        lino(8) = 3
        lino(9) = 13
        lino(10) = 20
        lino(11) = 16
        lino(12) = 12
        lino(13) = 9
        lino(14) = 17
        lino(15) = 19
        lino(16) = 11
        lino(17) = 14
        lino(18) = 18
        lino(19) = 15
        lino(20) = 10
        !write(*,*) '5'
! -------------------------------------------------------------------------------------------------
    elseif (conneo(5) .ne. 0 .and. conneo(6) .ne. 0 .and. &
            conneo(7) .ne. 0 .and. conneo(8) .ne. 0) then
!--------------------------------------------------------------------------------------------------
        lino(1) = 5
        lino(2) = 6
        lino(3) = 7
        lino(4) = 8
        lino(5) = 1
        lino(6) = 2
        lino(7) = 3
        lino(8) = 4
        lino(9) = 17
        lino(10) = 18
        lino(11) = 19
        lino(12) = 20
        lino(13) = 13
        lino(14) = 14
        lino(15) = 15
        lino(16) = 16
        lino(17) = 9
        lino(18) = 10
        lino(19) = 11
        lino(20) = 12
        !write(*,*) '6'
! -------------------------------------------------------------------------------------------------
    else
        ASSERT(.false.)
    end if
! -------------------------------------------------------------------------------------------------
    call jeecra(jexnum(conloc, indma), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(1)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(8)-1)
    zi(jconloc+5-1) = indno
    zi(jconloc+6-1) = zi(jmacsu+lino(12)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(13)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(20)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(16)-1)
    zi(jconloc+10-1) = indno+conneo(lino(4))
    zi(jconloc+11-1) = indno+conneo(lino(1))
    zi(jconloc+12-1) = indno+4+conneo(lino(1))
    zi(jconloc+13-1) = indno+4+conneo(lino(4))

    call jeecra(jexnum(conloc, indma+1), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma+1), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma+1), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(1)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+5-1) = indno
    zi(jconloc+6-1) = zi(jmacsu+lino(9)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(14)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(17)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(13)-1)
    zi(jconloc+10-1) = indno+conneo(lino(1))
    zi(jconloc+11-1) = indno+conneo(lino(2))
    zi(jconloc+12-1) = indno+4+conneo(lino(2))
    zi(jconloc+13-1) = indno+4+conneo(lino(1))

    call jeecra(jexnum(conloc, indma+2), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma+2), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma+2), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(2)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(3)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(7)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+5-1) = indno
    zi(jconloc+6-1) = zi(jmacsu+lino(10)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(15)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(18)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(14)-1)
    zi(jconloc+10-1) = indno+conneo(lino(2))
    zi(jconloc+11-1) = indno+conneo(lino(3))
    zi(jconloc+12-1) = indno+4+conneo(lino(3))
    zi(jconloc+13-1) = indno+4+conneo(lino(2))

    call jeecra(jexnum(conloc, indma+3), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma+3), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma+3), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(3)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(4)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(8)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(7)-1)
    zi(jconloc+5-1) = indno
    zi(jconloc+6-1) = zi(jmacsu+lino(11)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(16)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(19)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(15)-1)
    zi(jconloc+10-1) = indno+conneo(lino(3))
    zi(jconloc+11-1) = indno+conneo(lino(4))
    zi(jconloc+12-1) = indno+4+conneo(lino(4))
    zi(jconloc+13-1) = indno+4+conneo(lino(3))

    call jeecra(jexnum(conloc, indma+4), 'LONMAX', ival=13)
    call jeecra(jexnum(conloc, indma+4), 'LONUTI', ival=13)
    call jeveuo(jexnum(conloc, indma+4), 'E', jconloc)
    zi(jconloc+1-1) = zi(jmacsu+lino(5)-1)
    zi(jconloc+2-1) = zi(jmacsu+lino(6)-1)
    zi(jconloc+3-1) = zi(jmacsu+lino(7)-1)
    zi(jconloc+4-1) = zi(jmacsu+lino(8)-1)
    zi(jconloc+5-1) = indno
    zi(jconloc+6-1) = zi(jmacsu+lino(17)-1)
    zi(jconloc+7-1) = zi(jmacsu+lino(18)-1)
    zi(jconloc+8-1) = zi(jmacsu+lino(19)-1)
    zi(jconloc+9-1) = zi(jmacsu+lino(20)-1)
    zi(jconloc+10-1) = indno+4+conneo(lino(1))
    zi(jconloc+11-1) = indno+4+conneo(lino(2))
    zi(jconloc+12-1) = indno+4+conneo(lino(3))
    zi(jconloc+13-1) = indno+4+conneo(lino(4))

! -------------------------------------------------------------------------------------------------
    call jedema()
end subroutine
