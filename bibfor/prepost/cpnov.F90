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

subroutine cpnov(main, numa, coor, ind, nomnoe, conneo)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/reerel.h"
!
    integer(kind=8) ::  ind, numa, conneo(*)
    real(kind=8) :: coor(3, *)
    character(len=8) :: main
    character(len=24) :: nomnoe
!
!
! ----------------------------------------------------------------------
!         CREATION DU NOEUD SUPPLEMENTAIRE
!         SUR LE VOLUME 'DE LA ZONE DE CONTACT ESCLAVE'
! ----------------------------------------------------------------------
! IN        NUMA    NUMERO DE LA MAILLE COURANTE
! IN        IND     INDICE DU PREMIER NOEUD AJOUTE
! IN/JXVAR  NOMNOE  REPERTOIRE DE NOMS DES NOEUDS
! VAR       COOR    COORDONNEES DES NOEUDS
! IN        CONNEO  CONNECTIVIT2 DES ORDRE DES NOEUDS DU VOLUME
!                   ET DE LA FACE CORRESPONDANTE
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: lino(10), jtab, lgnd, iret
    integer(kind=8) :: inc1, inc2, aux
    real(kind=8) ::xe(3), xp(3), tabar(10*3)
!
    character(len=8) :: nomnd, eletyp
    character(len=24) :: valk
    character(len=16) :: knume

! ----------------------------------------------------------------------
!
! - INSERTION DU NOUVEAU NOEUD
!
!      NOM DU NOEUD CREE
    call codent(ind, 'G', knume)
    if (knume(1:1) == '*') then
        ASSERT(.false.)
    end if
    lgnd = lxlgut(knume)
    if (lgnd+1 .gt. 8) then
        call utmess('F', 'ALGELINE_16')
    end if
    nomnd = 'C'//knume(1:lgnd)
!
!      DECLARATION DU NOEUD CREE
    call jeexin(jexnom(nomnoe, nomnd), iret)
    if (iret .eq. 0) then
        call jecroc(jexnom(nomnoe, nomnd))
    else
        valk = nomnd
        call utmess('F', 'ALGELINE4_5', sk=valk)
    end if
!
! - CALCUL DES COORDONNEES DU NOUVEAU NOEUD
    call jeveuo(jexnum(main//'.CONNEX', numa), 'L', jtab)
    do inc1 = 1, 10
        lino(inc1) = zi(jtab+inc1-1)
    end do
    aux = 1
    do inc1 = 1, 10
        do inc2 = 1, 3
            tabar(aux+inc2-1) = coor(inc2, lino(inc1))
        end do
        aux = aux+3
    end do
    if (conneo(1) .eq. 0) then
        xe(1) = 1.d0/6.d0
        xe(2) = 1.d0/2.d0
        xe(3) = 1.d0/6.d0
    end if
    if (conneo(2) .eq. 0) then
        xe(1) = 1.d0/6.d0
        xe(2) = 1.d0/6.d0
        xe(3) = 1.d0/2.d0
    end if
    if (conneo(3) .eq. 0) then
        xe(1) = 1.d0/6.d0
        xe(2) = 1.d0/6.d0
        xe(3) = 1.d0/6.d0
    end if
    if (conneo(4) .eq. 0) then
        xe(1) = 1.d0/2.d0
        xe(2) = 1.d0/6.d0
        xe(3) = 1.d0/6.d0
    end if
    eletyp = 'T10'
    call reerel(eletyp, 10, 3, tabar, xe, xp)
    coor(1, ind) = xp(1)
    coor(2, ind) = xp(2)
    coor(3, ind) = xp(3)
!
end subroutine
