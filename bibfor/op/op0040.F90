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

subroutine op0040()
    implicit none
! OPERATEUR : INFO_RESU
! BUT       : FOURNIR LES COMPOSANTES DES CHAMPS PRESENTS DANS UNE SD
!             DE DONNEES RESULTAT
! ----------------------------------------------------------------------
! person_in_charge: nicolas.sellenet at edf.fr
#include "jeveux.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cmpcha.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ulexis.h"
#include "asterfort/ulopen.h"
!
    integer(kind=8) :: ifm, niv, ibid, nb_cmp, icmp, jatach, nbcham
    integer(kind=8) :: isy, iord, ifi, n2
    character(len=16) :: nomsym, nomfi
    character(len=19) :: resuin, nomcha
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
!
    call infniv(ifm, niv)
!
    ifi = 0
    nomfi = ' '
    call getvis(' ', 'UNITE', scal=ifi, nbret=n2)
    if (.not. ulexis(ifi)) then
        call ulopen(ifi, ' ', nomfi, 'NEW', 'O')
    end if
!
    call jemarq()
!
!     RECUPERATION DU NOM DU RESULTAT
    call getvid(' ', 'RESULTAT', scal=resuin, nbret=ibid)
!
    write (ifi, *) '-----------------------------------------------',&
     &                '------------'
    write (ifi, *)&
     & 'COMPOSANTES DES CHAMPS PRESENTS DANS LE RESULTAT : ', resuin
!
!     LECTURE DU NOMBRE DE CHAMPS PRESENTS ET DU NOMBRE D'ORDRE
    call jelira(resuin//'.DESC', 'NOMMAX', nbcham)
!
    do isy = 1, nbcham
        call jenuno(jexnum(resuin//'.DESC', isy), nomsym)
        call jenonu(jexnom(resuin//'.DESC', nomsym), ibid)
        call jeveuo(jexnum(resuin//'.TACH', ibid), 'L', jatach)
        iord = 1
        if (zk24(jatach-1+iord) (1:1) .ne. ' ') then
            write (ifi, *) '   - CHAMP ', nomsym, ' :'
            nomcha = zk24(jatach-1+iord) (1:19)
            call cmpcha(nomcha, cmp_name, cata_to_field, field_to_cata, nb_cmp)
            do icmp = 1, nb_cmp
                write (ifi, *) '      * ', cmp_name(icmp)
            end do
            AS_DEALLOCATE(vi=cata_to_field)
            AS_DEALLOCATE(vi=field_to_cata)
            AS_DEALLOCATE(vk8=cmp_name)
        end if
    end do
!
    write (ifi, *) '-----------------------------------------------',&
     &                '------------'
!
    call jedema()
end subroutine
