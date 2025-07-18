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

subroutine ss2mm2(mo, vecel, nomcas)
! INSPI  SS2MME
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
    character(len=8) :: mo, nomcas
    character(len=19) :: vecel
! ----------------------------------------------------------------------
!     BUT: TRAITER LE MOT-CLEF CAS_CHARGE DE LA COMMANDE
!          MACR_ELEM_STAT (POUR LES MACR_ELEM DU NIVEAU INFERIEUR)
!
!
!     IN:     MO : NOM DU MODELE
!          VECEL : NOM DU VECT_ELEM
!          NOMCAS: NOM DU CAS_CHARGE
!
!     OUT: VECEL EST  ENRICHI (EVENTUELLEMENT) DE L'OBJET .LISTE_CHAR
!
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: ma, nosma, nomacr
!
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ialsch
    integer(kind=8) :: iret, nbsma, nbssa
    character(len=24), pointer :: rerr(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
    call dismoi('NB_SS_ACTI', mo, 'MODELE', repi=nbssa)
    call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
!
    if (nbssa .eq. 0) goto 999
!
    call jeveuo(mo//'.MODELE    .SSSA', 'L', vi=sssa)
    call jeveuo(ma//'.NOMACR', 'L', vk8=vnomacr)
!
    call jeveuo(vecel//'.RERR', 'E', vk24=rerr)
    rerr(3) (1:3) = 'OUI'
!
    call jecrec(vecel//'.RELC', 'V V I', 'NO', 'CONTIG', 'CONSTANT', &
                1)
    call jeecra(vecel//'.RELC', 'LONMAX', nbsma)
    call jecroc(jexnom(vecel//'.RELC', nomcas))
    call jeveuo(jexnom(vecel//'.RELC', nomcas), 'E', ialsch)
!
!
!     -- REMPLISSAGE DE .RELC:
!     ------------------------------
!
!     -- ON VERIFIE QUE LES VECTEURS ELEMENTAIRES SONT CALCULES:
!     ----------------------------------------------------------
    do i = 1, nbsma
        if (sssa(i) .eq. 0) goto 3
        call jenuno(jexnum(ma//'.SUPMAIL', i), nosma)
        nomacr = vnomacr(i)
        call jeexin(jexnom(nomacr//'.LICA', nomcas), iret)
        if (iret .gt. 0) then
            zi(ialsch-1+i) = 1
        else
            zi(ialsch-1+i) = 0
        end if
3       continue
    end do
!
!
!
!
!
999 continue
    call jedema()
end subroutine
