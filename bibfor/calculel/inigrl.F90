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

subroutine inigrl(ligrel, igrel, nmax, adtabl, k24tab, &
                  nval)
! aslint: disable=
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
!
#include "asterfort/ini002.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
    character(len=*) :: ligrel
    integer(kind=8) :: igrel, nmax, adtabl(nmax), nval
    character(len=24) :: k24tab(nmax)
! ----------------------------------------------------------------------
!     BUT:
!     INITIALISER LE TYPE_ELEMENT ASSOCIE AU GREL  (INI00K)
!
!     IN:
!      LIGREL : NOM DU LIGREL A INITIALISER
!      IGREL  : NUMERO DU GREL
!      NMAX   : DIMENSION DE K24TAB ET ADTABL
!
!     OUT:
!         CREATION DES OBJETS JEVEUX PROPRES AU
!             TYPE_ELEMENT PRESENTS DANS LE LIGREL(IGREL).
!
!         NVAL  : NOMBRE DE NOMS RENDUS DANS K24TAB
!         K24TAB: TABLEAU DES NOMS DES OBJETS '&INEL.XXXX'
!         ADTABL : TABLEAU D'ADRESSES DES '&INEL.XXXXX'
!         ADTABL(I) = 0 SI L'OBJET CORRESPONDANT N'EXISTE PAS
!         SI NVAL > NMAX  : ON  S'ARRETE EN ERREUR FATALE
!
! ----------------------------------------------------------------------
    character(len=24) :: noliel
    character(len=16) :: nomte
    integer(kind=8) :: liel, l, te, k
!
!
!
    noliel = ligrel(1:19)//'.LIEL'
    call jeveuo(jexnum(noliel, igrel), 'L', liel)
    call jelira(jexnum(noliel, igrel), 'LONMAX', l)
    te = zi(liel-1+l)
    call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
!
!     -- ON MET LES ADRESSES A ZERO :
    do k = 1, nmax
        k24tab(k) = ' '
        adtabl(k) = 0
    end do
!
    call ini002(nomte, nmax, adtabl, k24tab, nval)
!
!
end subroutine
