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

subroutine foninf2(resu, typm, typfon, noma)
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    character(len=8) :: resu, typm, typfon, noma
!
!
!       ---------------------------------------------------------------
!       STOCKAGE D'INFOS UTILES DANS LA SD EN SORTIE DE DEFI_FOND_FISS
!       ---------------------------------------------------------------
!
! IN/OUT
!       RESU   : NOM DE LA SD EN SORTIE DE DEFI_FOND_FISS
!       TYPFON : TYPE DE FOND DE FISSURE (OUVERT ou FERME)
!       TYPM   : TYPE DE MAILLE DU FOND DE FISSURE (SEG2 ou SEG3)
!       NOMA   : NOM DU MAILLAGE
!
!
    integer(kind=8) ::  ibid, jinfo, n1, n2, n3
    character(len=8) :: syme, confin
    character(len=24) :: levsup, levinf
!
!     -----------------------------------------------------------------
!
    call jemarq()

    levsup = ''
    levinf = ''

    call jeveuo(noma//'.DIME', 'L', n1)
!
    call getvtx('LEVRE_SUP', 'GROUP_MA', iocc=1, nbval=0, nbret=n2)
    call getvtx('LEVRE_INF', 'GROUP_MA', iocc=1, nbval=0, nbret=n3)
!
    if (n2 .ne. 0) then
        call getvtx('LEVRE_SUP', 'GROUP_MA', iocc=1, scal=levsup)
    end if
!
    if (n3 .ne. 0) then
        call getvtx('LEVRE_INF', 'GROUP_MA', iocc=1, scal=levinf)
    end if
!
!     RECUPERATION DU MOT-CLE SYME
    call getvtx(' ', 'SYME', scal=syme, nbret=ibid)
    ASSERT(syme .eq. 'OUI' .or. syme .eq. 'NON')
!
!     RECUPERATION DU MOT-CLE CONFIG_INIT
    call getvtx(' ', 'CONFIG_INIT', scal=confin, nbret=ibid)
    ASSERT(confin .eq. 'DECOLLEE' .or. confin .eq. 'COLLEE')
!
!     CREATION DE L'OBJET .INFO DANS LA SD FOND_FISS
    call wkvect(resu//'.INFO', 'G V K24', 7, jinfo)
!
!     STOCKAGE DU MOT-CLE SYME
    zk24(jinfo-1+1) = syme
!
!     STOCKAGE DU MOT-CLE CONFIG_INIT
    zk24(jinfo-1+2) = confin
!
!     STOCKAGE DU MOT-CLE TYPE_FOND
    zk24(jinfo-1+3) = typfon
!
!     STOCKAGE DU NOM DU MAILLAGE
    zk24(jinfo-1+4) = noma

!     STOCKAGE DU TYPE DE MAILLE EN FOND DE FISSURE (3D)
    if (zi(n1-1+6) .eq. 3) then
        zk24(jinfo-1+5) = typm
    else
        zk24(jinfo-1+5) = ''
    end if

!     STOCKAGE DU NOM DE GROUPE DE MAIL LEVRESUP
    zk24(jinfo-1+6) = levsup

!     STOCKAGE DU NOM DE GROUPE DE MAIL LEVREINF
    zk24(jinfo-1+7) = levinf

    call jedema()
end subroutine
