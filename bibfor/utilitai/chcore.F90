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

subroutine chcore(chou)
    implicit none

!     BUT : TRANSFORMER UN CHAM_NO : REEL --> COMPLEXE

!     LE CHAMP COMPLEXE EST CONSTRUIT DE SORTE QUE:
!    - SA PARTIE REELLE CORRESPOND AUX VALEURS DU CHAMP REEL
!    - SA PARTIE IMAGINAIRE EST NULLE.
!     -----------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/getvid.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdchgd.h"
#include "asterfort/utmess.h"

    integer(kind=8) :: iret, jvale, nbval, jvalin, i, jrefe
    real(kind=8) :: zero
    parameter(zero=0.d0)
    character(len=3) :: tsca
    character(len=4) :: tych
    character(len=8) :: chou, chin, nomgd
    character(len=19) :: nume_equa, nume_equa_tmp, c19
    character(len=24) :: vale, valin, noojb
!---------------------------------------------------------------------

    call jemarq()

!   -- recuperation du champ reel
    call getvid(' ', 'CHAM_GD', scal=chin, nbret=iret)

!   -- verification : chin cham_no et reel ?
    call dismoi('TYPE_CHAMP', chin, 'CHAMP', repk=tych)
    if (tych .ne. 'NOEU') call utmess('F', 'UTILITAI_37', sk=chin)
    call dismoi('NOM_GD', chin, 'CHAMP', repk=nomgd)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    if (tsca .ne. 'R') then
        call utmess('F', 'UTILITAI_20', sk=chin)
    end if

!   -- copie chin --> chou
    call copisd('CHAMP', 'G', chin, chou)

    call dismoi("NUME_EQUA", chin, "CHAMP", repk=nume_equa)
    noojb = '12345678.NUMEC00000.PRNO'
    call gnomsd(chou, noojb, 15, 19)
    nume_equa_tmp = noojb(1:19)
    call copisd('NUME_EQUA', 'G', nume_equa, nume_equa_tmp)
    call jeveuo(nume_equa_tmp//'.REFN', 'E', jrefe)
    zk24(jrefe-1+2) = nomgd(1:5)//"C"
    c19 = chou
    call jeveuo(c19//'.REFE', 'E', jrefe)
    zk24(jrefe-1+2) = nume_equa_tmp

!   modifications de chou:
!   ======================

!   -- 1. ".vale"
!   --------------
    vale = chou
    vale(20:24) = '.VALE'

    call jelira(vale, 'LONMAX', nbval)
    call jedetr(vale)
    call jecreo(vale, 'G V C')
    call jeecra(vale, 'LONMAX', nbval)
    call jeveuo(vale, 'E', jvale)

    valin = vale
    valin(1:19) = chin
    call jeveuo(valin, 'L', jvalin)

    do i = 1, nbval
        zc(jvale+i-1) = dcmplx(zr(jvalin+i-1), zero)
    end do

!   -- 2. changement de la grandeur
!   ----------------------------
    call sdchgd(chou, 'C')

    call jedema()

end subroutine
