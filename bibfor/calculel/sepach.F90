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

subroutine sepach(carael, chinz, base, chreel, chimag)
!
    implicit none
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/cesvar.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedup1.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/nopar2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=1) :: base
    character(len=8) :: carael
    character(len=19) :: chreel, chimag
    character(len=*) :: chinz
!
!
!-----------------------------------------------------------------------
!
!
    integer(kind=8) :: gd, gdre, jdescr, jdesci, nbval, nbval2
    integer(kind=8) :: jvaler, jvalei, ivale, ier, iret1, iret2
    integer(kind=8) :: nmax1, nmax2, jncmpr, jncmpc, i, jrefe
    integer(kind=8) ::  nbsp
    character(len=8) :: nomgd, nomre
    character(len=4) :: typch
    character(len=1) :: ktyp
    character(len=19) :: canbva, chin, nume_equa, nume_equa_tmp
    character(len=24) :: ligrel, option, param, valk(2), noojb
    character(len=24), pointer :: celk(:) => null()
    complex(kind=8), pointer :: celv(:) => null()
    real(kind=8), pointer :: celvi(:) => null()
    real(kind=8), pointer :: celvr(:) => null()
    aster_logical :: l_complex
!
    call jemarq()
!
    chin = chinz
!
    call dismoi("DOCU", chin, "CHAMP", repk=typch)
    call dismoi("NUM_GD", chin, "CHAMP", repi=gd)
!
    ier = 0
!
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nomgd)
    if ((nomgd(7:7) .ne. ' ') .or. (nomgd(5:6) .ne. '_C')) then
        call utmess('F', 'CALCULEL4_80', sk=nomgd)
    end if
    nomre = nomgd(1:4)//'_R'
    call jenonu(jexnom('&CATA.GD.NOMGD', nomre), gdre)
!
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nmax1)
    call jelira(jexnum('&CATA.GD.NOMCMP', gdre), 'LONMAX', nmax2)
!
    if (nmax1 .ne. nmax2) then
        valk(1) = nomgd
        valk(2) = nomre
        call utmess('F', 'CALCULEL4_81', nk=2, valk=valk)
    end if
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gdre), 'L', jncmpr)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', jncmpc)
!
    do i = 1, nmax1
        if (zk8(jncmpr-1+i) .ne. zk8(jncmpc-1+i)) then
            ier = 1
            goto 10
        end if
10      continue
    end do
!
    if (ier .ne. 0) then
        valk(1) = nomgd
        valk(2) = nomre
        call utmess('F', 'CALCULEL4_82', nk=2, valk=valk)
    end if
!
!     -- CHAM_NO :
!     -------------------
    if (typch .eq. 'CHNO') then
        noojb = '12345678.NUMEC00000.PRNO'
        call gnomsd(chreel, noojb, 15, 19)
        nume_equa_tmp = noojb(1:19)
        call dismoi("NUME_EQUA", chin, "CHAMP", repk=nume_equa)
        call copisd('NUME_EQUA', 'G', nume_equa, nume_equa_tmp)
        call jeveuo(nume_equa_tmp//'.REFN', 'E', jrefe)
        zk24(jrefe-1+2) = nomre

        call jedup1(chin//'.REFE', base, chreel//'.REFE')
        call jeveuo(chreel//'.REFE', 'E', jrefe)
        zk24(jrefe-1+2) = nume_equa_tmp
!
        call jedup1(chreel//'.REFE', base, chimag//'.REFE')
!
        call jelira(chin//'.VALE', 'LONMAX', nbval)
        call jeveuo(chin//'.VALE', 'L', ivale)
        call jelira(chin//'.VALE', 'TYPE', cval=ktyp)
        l_complex = (ktyp == 'C')
!
        call wkvect(chreel//'.VALE', base//' V R', nbval, jvaler)
        call wkvect(chimag//'.VALE', base//' V R', nbval, jvalei)
!
        if (l_complex) then
            do i = 1, nbval
                zr(jvaler-1+i) = dble(zc(ivale-1+i))
                zr(jvalei-1+i) = dimag(zc(ivale-1+i))
            end do
        else
            do i = 1, nbval
                zr(jvaler-1+i) = zr(ivale-1+i)
                zr(jvalei-1+i) = 0.D0
            end do
        end if
!
!
!
!     -- CHAM_ELEM :
!     -------------------
    else if (typch .eq. 'CHML') then
        call jeveuo(chin//'.CELK', 'L', vk24=celk)
        ligrel = celk(1)
        option = celk(2)
!
        call nopar2(option, nomre, 'OUT', param)
!
!       -- SI LE CHIN A DES SOUS-POINTS, IL FAUT ALLOUER CHREEL
!          ET CHIMAG AVEC DES SOUS-POINTS :
        call dismoi('MXNBSP', chin, 'CHAM_ELEM', repi=nbsp)
        if (nbsp .gt. 1) then
            canbva = '&&SEPACH.CANBVA'
            call cesvar(carael, ' ', ligrel, canbva)
            call alchml(ligrel, option, param, base, chreel, &
                        iret1, canbva)
            call alchml(ligrel, option, param, base, chimag, &
                        iret2, canbva)
            call detrsd('CHAM_ELEM_S', canbva)
        else
            call alchml(ligrel, option, param, base, chreel, &
                        iret1, ' ')
            call alchml(ligrel, option, param, base, chimag, &
                        iret2, ' ')
        end if
!
        ASSERT((iret1 .eq. 0) .or. (iret2 .eq. 0))
!
        call jelira(chin//'.CELV', 'LONMAX', nbval)
        call jeveuo(chin//'.CELV', 'L', vc=celv)
!
        call jelira(chreel//'.CELV', 'LONMAX', nbval2)
        if (nbval2 .ne. nbval) then
            ASSERT(nbval .gt. nbval2)
            call juveca(chreel//'.CELV', nbval)
            call juveca(chimag//'.CELV', nbval)
        end if
!
!
        call jeveuo(chreel//'.CELV', 'E', vr=celvr)
        call jeveuo(chimag//'.CELV', 'E', vr=celvi)
!
!
        do i = 1, nbval
            celvr(i) = dble(celv(i))
            celvi(i) = dimag(celv(i))
        end do
!
!     -- CART :
!     -------------------
    else if (typch .eq. 'CART') then
        call jedup1(chin//'.DESC', base, chreel//'.DESC')
        call jedup1(chin//'.NOMA', base, chreel//'.NOMA')
        call jedup1(chin//'.NOLI', base, chreel//'.NOLI')
        call jedup1(chin//'.LIMA', base, chreel//'.LIMA')
        call jeveuo(chreel//'.DESC', 'E', jdescr)
        zi(jdescr) = gdre
!
        call jedup1(chin//'.DESC', base, chimag//'.DESC')
        call jedup1(chin//'.NOMA', base, chimag//'.NOMA')
        call jedup1(chin//'.NOLI', base, chimag//'.NOLI')
        call jedup1(chin//'.LIMA', base, chimag//'.LIMA')
        call jeveuo(chimag//'.DESC', 'E', jdesci)
        zi(jdesci) = gdre
!
        call jelira(chin//'.VALE', 'LONMAX', nbval)
        call jeveuo(chin//'.VALE', 'L', ivale)
!
        call wkvect(chreel//'.VALE', base//' V R', nbval, jvaler)
        call wkvect(chimag//'.VALE', base//' V R', nbval, jvalei)
!
        do i = 1, nbval
            zr(jvaler-1+i) = dble(zc(ivale-1+i))
            zr(jvalei-1+i) = dimag(zc(ivale-1+i))
        end do
!
    end if
!
    call jedema()
!
end subroutine
