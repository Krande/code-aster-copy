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

subroutine chpeva(chou)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
!     BUT : TRAITER :
!      -  OPTION:'EVAL' DE LA COMMANDE CREA_CHAMP
!      -  COMMANDE CALC_CHAMP
!     ATTENTION: CETTE ROUTINE N'EST PAS UN UTILITAIRE :
!                ELLE FAIT DES CALL GETVXX.
!     -----------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/ceseva.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnseva.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: n1, ib, npara, k, nncp, ibid
    character(len=4) :: typ1, typ2, knum
    character(len=8) :: chou, nomgd
    character(len=19) :: ligrel, chs1, chs2, chins, chin
    character(len=24), pointer :: lpara1(:) => null()
    character(len=24), pointer :: lpara2(:) => null()
!     -----------------------------------------------------------------
!
    call jemarq()
!
!
    chins = '&&CHPEVA.CHINS'
    chs2 = '&&CHPEVA.CHS2'
!
!
! 1. CALCUL DE:
!      CHIN  : CHAMP A EVALUER (DE FONCTIONS)
!      NOMGD : GRANDEUR ASSOCIEE A CHIN
!      .LPARA1: NOMS DES CHAMPS PARAMETRES POUR LES FONCTIONS
! ------------------------------------------------------------------
!
    call getvid(' ', 'CHAM_F', scal=chin, nbret=ib)
!
    call dismoi('NOM_GD', chin, 'CHAMP', repk=nomgd)
    if (nomgd .ne. 'NEUT_F') then
        call utmess('F', 'MODELISA4_13')
    end if
!
    call getvid(' ', 'CHAM_PARA', nbval=0, nbret=n1)
    npara = -n1
    AS_ALLOCATE(vk24=lpara1, size=npara)
    call getvid(' ', 'CHAM_PARA', nbval=npara, vect=lpara1, nbret=n1)
!
!
! 2.  ON VERIFIE QUE LES CHAMPS PARA ONT LA MEME DISCRETISATION:
!     ET ON LES TRANSFORME EN CHAMPS SIMPLES
!     CALCUL DE .LPARA2: NOMS DES CHAMP_S PARAMETRES POUR LES FONCTIONS
! ------------------------------------------------------------
    AS_ALLOCATE(vk24=lpara2, size=npara)
    call dismoi('TYPE_CHAMP', chin, 'CHAMP', repk=typ1)
    do k = 1, npara
        call dismoi('TYPE_CHAMP', lpara1(k), 'CHAMP', repk=typ2)
        if (typ1 .ne. typ2) then
            call utmess('F', 'MODELISA4_14')
        end if
!
        call codent(k, 'G', knum)
        chs1 = '&&CHPEVA.'//knum
        if (typ1 .eq. 'NOEU') then
            call cnocns(lpara1(k), 'V', chs1)
!
        else if (typ1 .eq. 'CART') then
            call carces(lpara1(k), 'ELEM', ' ', 'V', chs1, &
                        'A', ib)
!
        else if (typ1(1:2) .eq. 'EL') then
            call celces(lpara1(k), 'V', chs1)
        end if
        lpara2(k) = chs1
    end do
!
!
! 3.  -- ON APPELLE LA ROUTINE D'EVAL APPROPRIEE :
! ------------------------------------------------------------
    ASSERT((typ1 .eq. 'NOEU') .or. (typ1(1:2) .eq. 'EL'))
    if (typ1 .eq. 'NOEU') then
        call cnocns(chin, 'V', chins)
        call cnseva(chins, npara, lpara2, chs2)
        call cnscno(chs2, ' ', 'NON', 'G', chou, &
                    'F', ibid)
        call detrsd('CHAM_NO_S', chins)
        call detrsd('CHAM_NO_S', chs2)
!
    else if (typ1(1:2) .eq. 'EL') then
        call celces(chin, 'V', chins)
        call ceseva(chins, npara, lpara2, chs2)
        call dismoi('NOM_LIGREL', chin, 'CHAMP', repk=ligrel)
        call cescel(chs2, ligrel, ' ', ' ', 'NON', &
                    nncp, 'G', chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', chins)
        call detrsd('CHAM_ELEM_S', chs2)
!
    end if
!
!
! 7. MENAGE :
! -----------------------------------------------------
    AS_DEALLOCATE(vk24=lpara1)
    do k = 1, npara
        if (typ1 .eq. 'NOEU') then
            call detrsd('CHAM_NO_S', lpara2(k))
!
        else
            call detrsd('CHAM_ELEM_S', lpara2(k))
        end if
    end do
    AS_DEALLOCATE(vk24=lpara2)
!
    call jedema()
!
end subroutine
