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

subroutine dbg_base()
    implicit none
!
! person_in_charge: mathieu.courtois@edf.fr
!
! Check that the splitting of the results databases works as expected.
!
! This subroutine is only used by unittest and enabled using:
!   DEBUT(DEBUG=_F(VERI_BASE='GLOBALE' or 'VOLATILE'))
!
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelibe.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: exist
    integer(kind=8), parameter :: n = 5
    integer(kind=8) :: lfic, mfic
    common/fenvje/lfic(n), mfic
!
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
!
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &               kitlec(n), kitecr(n), kiadm(n),&
     &               iitlec(n), iitecr(n), nitecr(n), kmarq(n)
!
    character(len=1) :: dblett
    character(len=4) :: dbname, valk(2)
    character(len=16) :: answer
    character(len=19) :: objname
    integer(kind=8) :: dbsize, recsize, recint
    integer(kind=8) :: nbrec, objint, nbobj, ic, i, total, perc
    integer(kind=8) :: jadr, vali(10), more
!
    answer = ' '
    call getvtx('DEBUG', 'VERI_BASE', iocc=1, scal=answer)
    call getvis('DEBUG', 'VERI_BASE_NB', iocc=1, scal=perc)
    dblett = answer(1:1)
!
    ic = index(classe, dblett)
    ASSERT(ic .ne. 0)
    dbname = nomfic(ic)
    recint = longbl(ic)*1024
    recsize = 8*recint
    nbrec = nblmax(ic)
    dbsize = recsize*nbrec
    objint = recint/4
!   with 2 TB, noobj=10486101
    nbobj = 4*nbrec*perc/100
    vali(1) = dbsize
    vali(10) = mfic*1024
    vali(9) = lfic(ic)*1000
    vali(2) = recsize
    vali(3) = recint
    vali(4) = nbrec
    vali(5) = objint
    vali(6) = nbobj
    vali(7) = nbobj*objint*8
    vali(8) = nbobj*objint
    call utmess('I', 'JEVEUX1_97', vali=vali, ni=10)

    inquire (file=dbname//'.1', exist=exist)
    ASSERT(exist)
    inquire (file=dbname//'.2', exist=exist)
    ASSERT(.not. exist)

    i = 0
    total = 0
    more = 10
    do while ((i < nbobj .and. .not. exist) .or. (exist .and. more .gt. 0))
        i = i+1
        write (objname, '(A12,I7)') '&&VERI_BASE.', i
        total = total+objint
        call wkvect(objname, dblett//' V I', objint, jadr)
        call jelibe(objname)
        inquire (file=dbname//'.2', exist=exist)
        if (exist) then
            more = more-1
        end if
        call utmess('I', 'JEVEUX1_98', sk=objname, si=total)
    end do

!   TEST_RESU
    inquire (file=dbname//'.2', exist=exist)
    valk(1) = dbname
    if (exist) then
        valk(2) = "OK"
    else
        valk(2) = "NOOK"
    end if
    call utmess('I', 'JEVEUX1_99', nk=2, valk=valk)
end subroutine
