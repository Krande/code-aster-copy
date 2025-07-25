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

subroutine cclpco(option, resuou, numord, nbpaou, lipaou, &
                  lichou)
    implicit none
!     --- ARGUMENTS ---
#include "jeveux.h"
!
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/assert.h"
    integer(kind=8) :: nbpaou, numord
    character(len=8) :: resuou
    character(len=8) :: lipaou(*)
    character(len=16) :: option
    character(len=24) :: lichou(*)
!  CALC_CHAMP - DETERMINATION LISTE DE PARAMETRES ET LISTE DE CHAMPS OUT
!  -    -                     -        -                      -      -
! ----------------------------------------------------------------------
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION A CALCULER
!   RESUOU  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT OUT
!   NUMORD  I    NUMERO D'ORDRE COURANT
!
! OUT :
!   NBPAOU  I    NOMBRE DE PARAMETRES OUT
!   LIPAOU  K8*  LISTE DES PARAMETRES OUT
!   LICHOU  K8*  LISTE DES CHAMPS OUT
! ----------------------------------------------------------------------
! person_in_charge: nicolas.sellenet at edf.fr
    integer(kind=8) :: opt, iaopds, iapara, nparin, ipara, ierd
    integer(kind=8) :: nparou, kpara, nugd
!
    character(len=8) :: nomgd, tsca
    character(len=19) :: nochou
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), opt)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iaopds)
    call jeveuo(jexnum('&CATA.OP.OPTPARA', opt), 'L', iapara)
!
    nparin = zi(iaopds-1+2)
    nparou = zi(iaopds-1+3)
!
    if (nparou .eq. 1) then
        ipara = 1

    elseif (nparou .eq. 2) then
!       -- on cherche le parametre de type reel :
        ipara = 0
        do kpara = 1, 2
            nugd = zi(iaopds-1+4+nparin+kpara)
            call jenuno(jexnum('&CATA.GD.NOMGD', nugd), nomgd)
            call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
            if (tsca .eq. 'R') ipara = kpara
        end do
        ASSERT(ipara .gt. 0)
    else
        ASSERT(.false.)
    end if

    nbpaou = 1
    lipaou(1) = zk8(iapara+nparin+ipara-1)
    call rsexch(' ', resuou, option, numord, nochou, ierd)
    lichou(1) = nochou

    call jedema()

end subroutine
