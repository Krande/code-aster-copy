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

subroutine agcart(ngdmxn, chinz)
    implicit none
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/nbec.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ngdmxn
    character(len=19) :: chin
    character(len=*) :: chinz
! --------------------------------------------------------------------
! Agrandissement de la carte chin, ngdmxn etant le nouveau nombre
!   maximum de couples (entite,valeur) a stocker
!
!   Attention : agcart n'agrandit pas l'objet .lima
!               celui-ci est agrandi directement par nocart
!     il est donc dangereux d'appeler agcart en dehors de nocart
! --------------------------------------------------------------------
!  ngdmxn       - in     - i    - : nouveau nombre max de couples
!               -        -      -   (entite,valeur) a stocker
! --------------------------------------------------------------------
!  chinz        - in     - k*(*)- : nom de la carte a redimensionner -
!               - jxvar  -      -   on realloue et on recopie leurs
!               -        -      -   anciennes valeurs pour les objets-
!               -        -      -   chin.desc
!               -        -      -   chin.vale
!               -        -      -   chin.noma
!               -        -      -   chin.noli
! --------------------------------------------------------------------
    character(len=1) :: base
    character(len=24) :: descav
    integer(kind=8) :: jdesca, jdesc, nec, iec, ngdmxa, nedit, ied, ideca, idec
    integer(kind=8) ::  ncmp, igd
! ----------------------------------------------------------------------
    call jemarq()
    chin = chinz
    call jelira(chin//'.DESC', 'CLAS', cval=base)
!
!
! --- AGRANDISSEMENT DE .DESC:
! ------------------------------
    descav = '&&AGCART.DESCAV'
    call jedupo(chin//'.DESC', 'V', descav, .false._1)
    call jeveuo(descav, 'E', jdesca)
    igd = zi(jdesca-1+1)
    nec = nbec(igd)
    ngdmxa = zi(jdesca-1+2)
    nedit = zi(jdesca-1+3)
    ASSERT(ngdmxn .gt. ngdmxa)
!
    call jedetr(chin//'.DESC')
    call wkvect(chin//'.DESC', base//' V I', 3+ngdmxn*(2+nec), jdesc)
    call jeecra(chin//'.DESC', 'DOCU', cval='CART')
!
    zi(jdesc-1+1) = igd
    zi(jdesc-1+2) = ngdmxn
    zi(jdesc-1+3) = nedit
!
    do ied = 1, nedit
        zi(jdesc-1+3+(ied-1)*2+1) = zi(jdesca-1+3+(ied-1)*2+1)
        zi(jdesc-1+3+(ied-1)*2+2) = zi(jdesca-1+3+(ied-1)*2+2)
    end do
!
    do ied = 1, nedit
        ideca = 3+2*ngdmxa+nec*(ied-1)
        idec = 3+2*ngdmxn+nec*(ied-1)
        do iec = 1, nec
            zi(jdesc-1+idec+iec) = zi(jdesca-1+ideca+iec)
        end do
    end do
    call jedetr(descav)
!
!
!
! ---  AGRANDISSEMENT DE VALE:
! ------------------------------
    call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmp)
    call juveca(chin//'.VALE', ngdmxn*ncmp)
    call jeecra(chin//'.VALE', 'LONUTI', ngdmxn*ncmp)
!
!
! ---  AGRANDISSEMENT DE NOLI
! ------------------------------
    call juveca(chin//'.NOLI', ngdmxn)
    call jeecra(chin//'.NOLI', 'LONUTI', ngdmxn)
!
!
! ---  AGRANDISSEMENT DE LIMA : ON NE FAIT RIEN :
!      C'EST NOCART QUI AGRANDIT .LIMA SI NECESSAIRE
!
!
    call jedema()
!
end subroutine
