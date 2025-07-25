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

subroutine trcart(ific, nocc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utcmp1.h"
#include "asterfort/utmess.h"
#include "asterfort/tresu_carte.h"
#include "asterfort/tresu_ordgrd.h"
#include "asterfort/tresu_read_refe.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/tresu_tole.h"
#include "asterfort/int_to_char8.h"
    integer(kind=8), intent(in) :: ific
    integer(kind=8), intent(in) :: nocc
!     COMMANDE:  TEST_RESU
!                MOT CLE FACTEUR "CARTE"
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: iocc, iret
    integer(kind=8) :: n1, n2, n3, n1r, n2r, n3r, ivari, n1a, n1b
    integer(kind=8) :: nl1, nl2, vali, valir
    integer(kind=8) :: jnuma
    real(kind=8) :: epsi, epsir, valr, valrr
    complex(kind=8) :: valc, valcr
    character(len=1) :: typres
    character(len=16) :: nom_vari
    character(len=8) :: crit, noddl, nomma, nomail, nomgd
    character(len=11) :: motcle
    character(len=19) :: cham19
    character(len=16) :: tbtxt(2), tbref(2)
    character(len=24) :: nogrma
    character(len=200) :: lign1, lign2
    aster_logical :: lref, l_parallel_mesh
    aster_logical :: skip
    real(kind=8) :: ordgrd
!     ------------------------------------------------------------------
    call jemarq()
!
    motcle = 'CARTE'
!
    do iocc = 1, nocc
        lign1 = ' '
        lign2 = ' '
        noddl = ' '
        call getvid('CARTE', 'CHAM_GD', iocc=iocc, scal=cham19, nbret=n1)
        lign1(1:21) = '---- '//motcle(1:9)
        lign1(22:22) = '.'
        lign2(1:21) = '     '//cham19(1:8)
        lign2(22:22) = '.'
        call tresu_read_refe('CARTE', iocc, tbtxt)
!
        call getvtx('CARTE', 'NOM_CMP', iocc=iocc, scal=noddl, nbret=n1)
        ASSERT(n1 .eq. 1)
        nl1 = lxlgut(lign1)
        nl2 = lxlgut(lign2)
        lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_CMP'
        lign2(1:nl2+16) = lign2(1:nl2-1)//' '//noddl
        lign1(nl1+17:nl1+17) = '.'
        lign2(nl2+17:nl2+17) = '.'
!
!
        call tresu_tole(epsi, mcf='CARTE', iocc=iocc)
        call getvtx('CARTE', 'CRITERE', iocc=iocc, scal=crit, nbret=n1)
!
        call getvr8('CARTE', 'VALE_CALC', iocc=iocc, scal=valr, nbret=n1)
        call getvis('CARTE', 'VALE_CALC_I', iocc=iocc, scal=vali, nbret=n2)
        call getvc8('CARTE', 'VALE_CALC_C', iocc=iocc, scal=valc, nbret=n3)
!
        skip = .false.
        ordgrd = 1.d0
        if (n1 .eq. 1) then
            typres = 'R'
            call tresu_ordgrd(valr, skip, ordgrd, mcf='CARTE', iocc=iocc)
        else if (n2 .eq. 1) then
            typres = 'I'
        else
            ASSERT(n3 .eq. 1)
            typres = 'C'
        end if
! ----------------------------------------------------------------------
        lref = .false.
        call getvr8('CARTE', 'PRECISION', iocc=iocc, scal=epsir, nbret=iret)
        if (iret .ne. 0) then
            lref = .true.
            call getvr8('CARTE', 'VALE_REFE', iocc=iocc, scal=valrr, nbret=n1r)
            call getvis('CARTE', 'VALE_REFE_I', iocc=iocc, scal=valir, nbret=n2r)
            call getvc8('CARTE', 'VALE_REFE_C', iocc=iocc, scal=valcr, nbret=n3r)
            ASSERT(n1r .eq. n1 .and. n2r .eq. n2 .and. n3r .eq. n3)
        end if
        if (skip .and. .not. lref) then
            call utmess('I', 'TEST0_11')
        end if
! ----------------------------------------------------------------------
!
!
        call getvtx('CARTE', 'NOM_CMP', iocc=iocc, scal=noddl, nbret=n1)
        call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=nomma)
        l_parallel_mesh = isParallelMesh(nomma)
        if (l_parallel_mesh) then
            call utmess('F', 'MODELISA7_88')
        end if
!
        call getvem(nomma, 'MAILLE', 'CARTE', 'MAILLE', iocc, &
                    1, nomail, n1a)
        if (n1a .eq. 0) then
            call getvem(nomma, 'GROUP_MA', 'CARTE', 'GROUP_MA', iocc, &
                        1, nogrma, n1b)
            ASSERT(n1b .eq. 1)
            call jelira(jexnom(nomma//'.GROUPEMA', nogrma), 'LONUTI', ival=n1b)
            if (n1b .ne. 1) call utmess('F', 'TEST0_20', sk=nogrma, si=n1b)
            call jeveuo(jexnom(nomma//'.GROUPEMA', nogrma), 'L', jnuma)
            nomail = int_to_char8(zi(jnuma))
        else
            if (l_parallel_mesh) then
                call utmess('F', 'MODELISA7_86')
            end if
        end if
!
        if (n1a .ne. 0) then
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' MAILLE'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nomail
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
        else
            nl1 = lxlgut(lign1)
            nl2 = lxlgut(lign2)
            lign1(1:nl1+16) = lign1(1:nl1-1)//' GROUP_MA'
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//nogrma
            lign1(nl1+17:nl1+17) = '.'
            lign2(nl2+17:nl2+17) = '.'
        end if
!
!
        call dismoi('NOM_GD', cham19, 'CHAMP', repk=nomgd)
        call utcmp1(nomgd, 'CARTE', iocc, noddl, ivari, nom_vari)
!
        nl1 = lxlgut(lign1)
        nl1 = lxlgut(lign1(1:nl1-1))
        nl2 = lxlgut(lign2)
        nl2 = lxlgut(lign2(1:nl2-1))
        write (ific, *) lign1(1:nl1)
        write (ific, *) lign2(1:nl2)
!
        if (lref) then
            tbref(1) = tbtxt(1)
            tbref(2) = tbtxt(2)
            tbtxt(1) = 'NON_REGRESSION'
        end if
        call tresu_carte(cham19, nomail, noddl, tbtxt, vali, &
                         valr, valc, typres, epsi, crit, &
                         .true._1, ignore=skip, compare=ordgrd)
        if (lref) then
            call tresu_carte(cham19, nomail, noddl, tbref, valir, &
                             valrr, valcr, typres, epsir, crit, &
                             .false._1)
        end if
        write (ific, *) ' '
!
    end do
!
!
    call jedema()
end subroutine
