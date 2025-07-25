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

subroutine u195tb(chou)
!
!
!     TRAITEMENT DE COMMANDE:   CREA_CHAMP / OPTION: 'EXTR' / TABLE
!
!     " CREATION D'UN CHAMP A PARTIR D'UNE TABLE "
!
    implicit none
!
!     ------------------------------------------------------------------
! 0.1. ==> ARGUMENT
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cescar.h"
#include "asterfort/cescel.h"
#include "asterfort/cnscno.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tabchs.h"
#include "asterfort/utmess.h"
    character(len=*) :: chou
!
! 0.2. ==> COMMUNS
!
!
!      ==> VARIABLES LOCALES
!
    integer(kind=8) :: n1, n2, ibid, nncp
    character(len=3) :: prol0
    character(len=8) :: nomgd, ma, mo
    character(len=16) :: tychlu, option, typchs, typch2
    character(len=19) :: chs, tabin, ligrel
    character(len=8), pointer :: lgrf(:) => null()
!
    call jemarq()
!
    chs = '&&U195TB.CHAMP_S'
!
    call getvid(' ', 'TABLE', scal=tabin, nbret=n1)
    call getvtx(' ', 'TYPE_CHAM', scal=tychlu, nbret=n1)
!
    typchs = tychlu(1:4)
    typch2 = typchs
    if (typch2 .eq. 'CART') typch2 = 'ELEM'
    nomgd = tychlu(6:11)
!
!     VERIFICATIONS
    if (typchs .eq. 'NOEU' .or. typchs .eq. 'CART') then
        call getvid(' ', 'MAILLAGE', scal=ma, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'MODELISA7_61')
        end if
        option = ' '
        mo = ' '
    else if (typchs(1:2) .eq. 'EL') then
        call getvid(' ', 'MODELE', nbval=0, nbret=n1)
        call getvtx(' ', 'OPTION', nbval=0, nbret=n2)
        if (n1 .ne. 0) then
            n1 = -n1
            call getvid(' ', 'MODELE', scal=mo, nbret=n1)
        end if
        if (n2 .ne. 0) then
            n2 = -n2
            call getvtx(' ', 'OPTION', scal=option, nbret=n2)
        end if
        if (n1 .eq. 0 .or. n2 .eq. 0) then
            call utmess('F', 'MODELISA7_62')
        end if
        call jeveuo(mo//'.MODELE    .LGRF', 'L', vk8=lgrf)
        ma = lgrf(1)
    end if
!
!     CREATION DU CHAMP SIMPLE
    call tabchs(tabin, typch2, 'V', nomgd, ma, &
                chs)
!
!     TRANSFORMATION : CHAM_S --> CHAMP
    call getvtx(' ', 'PROL_ZERO', scal=prol0, nbret=n1)
    if (n1 .eq. 0) prol0 = 'NON'
    if (typchs .eq. 'NOEU') then
        call cnscno(chs, ' ', prol0, 'G', chou, &
                    'F', ibid)
    else if (typchs(1:2) .eq. 'EL') then
        call dismoi('NOM_LIGREL', mo, 'MODELE', repk=ligrel)
!        -- POUR L'UTILISATEUR DISTRAIT :
        if (option .eq. 'VARI_ELGA') option = 'RAPH_MECA'
        call cescel(chs, ligrel, option, ' ', prol0, &
                    nncp, 'G', chou, 'F', ibid)
    else if (typchs .eq. 'CART') then
        call cescar(chs, chou, 'G')
    else
        ASSERT(.false.)
    end if
    call jedema()
!
end subroutine
