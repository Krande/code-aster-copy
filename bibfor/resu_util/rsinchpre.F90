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

subroutine rsinchpre(nomsd, nomch, acces, ier)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxliis.h"
!
    character(len=*), intent(in) :: nomsd, nomch, acces
    integer(kind=8), intent(out) :: ier
!      VERIFICATION DE FAISABILITE
!      INTERPOLATION D'UN CHAMP_19 A PARTIR D'1 SD RESULTAT-COMPOSE
! ----------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA STRUCTURE "RESULTAT"
! IN  : NOMCH  : NOM SYMBOLIQUE DU CHAMP CHERCHE.
! IN  : ACCES  : NOM SYMBOLIQUE DE LA VARIABLE D'ACCES.
!
! OUT : IER    : CODE_RETOUR :
!                10 --> IL N'EXISTE AUCUN CHAMP POUR L'INTERPOLATION.
!                20 --> LA VARIABLE D'ACCES EST ILLICITE.
!                30 --> LE CHAMP A INTERPOLER EST ILLICITE.
! ----------------------------------------------------------------------
    character(len=4) :: type, tysca
    character(len=8) :: nomobj, k8debu, k8maxi, k8ent
    character(len=19) :: noms2
    character(len=16) :: acce2, nomc2
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacces, iatach
    integer(kind=8) :: iatava, iadesc, idebu, ier1, ier2
    integer(kind=8) :: iloty, imaxi, nbordr
    aster_logical :: vide

!-----------------------------------------------------------------------
    call jemarq()
    acce2 = acces
    noms2 = nomsd
    nomc2 = nomch
    ier = 0
!
!     -- VERIFICATION DE LA VARIABLE D'ACCES:
!     ---------------------------------------
!

    if (ier .eq. 0) then
        call jenonu(jexnom(noms2//'.NOVA', acce2), iacces)
        if (iacces .eq. 0) then
            ier = 20
        end if
    end if

    if ((ier .eq. 0) .and. (iacces .ne. 0)) then
        call jeveuo(jexnum(noms2//'.TAVA', iacces), 'L', iatava)
        nomobj = zk8(iatava-1+1)
        k8debu = zk8(iatava-1+2)
        call lxliis(k8debu, idebu, ier1)
        k8maxi = zk8(iatava-1+3)
        call lxliis(k8maxi, imaxi, ier2)
        if ((ier1 .ne. 0) .or. (ier2 .ne. 0)) then
            ier = 20
        else
            call jelira(noms2//nomobj, 'TYPE', cval=type)
            call jelira(noms2//nomobj, 'LTYP', iloty)
            call codent(iloty, 'G', k8ent)
            tysca = type(1:1)//k8ent(1:3)
            if (tysca .ne. 'R8  ') then
                ier = 20
            end if
        end if
    end if
!
!     -- VERIFICATION DU NOM DE CHAMP:
!     --------------------------------
!
    if (ier .eq. 0) then
        call jenonu(jexnom(noms2//'.DESC', nomc2), iadesc)
        if (iadesc .eq. 0) then
            ier = 21
        end if
    end if
!
!     -- VERIFICATION DES CHAMPS EXISTANTS
!     ------------------------------------
!
    if ((ier .eq. 0) .and. (iacces .ne. 0) .and. (iadesc .ne. 0)) then
        call jelira(noms2//'.ORDR', 'LONUTI', nbordr)
        call jeveuo(jexnum(noms2//'.TACH', iadesc), 'L', iatach)
        vide = .true.
        do i = 1, nbordr
            if (zk24(iatach-1+i) (1:1) .ne. ' ') then
                vide = .false.
            end if
        end do
        if (vide) then
            ier = 10
        end if
    end if
!
    call jedema()
end subroutine rsinchpre
