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
!
subroutine rfmge1(modgen)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: modgen
!
!     OPERATEUR "RECU_FONCTION"  MOT CLE "RESU_GENE"
!                                CONCEPT MODE_GENE
!     ------------------------------------------------------------------
    integer(kind=8) :: n1, ncmp, iret, jordr, lpro, lvar, lfon, nbordr, im, iord, iad
    integer(kind=8) :: nbmode, i, istru
    real(kind=8) :: epsi
    character(len=4) :: interp(2)
    character(len=8) :: k8b, crit, mode
    character(len=14) :: nugene
    character(len=16) :: k16b, nomcmd, typcon, nomcha, npara
    character(len=19) :: noch19, nomfon, knume
    integer(kind=8), pointer :: nequ(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    character(len=24), pointer :: refe(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call getres(nomfon, typcon, nomcmd)
!
    interp(1) = 'LIN '
    interp(2) = 'LIN '
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=n1)
    if (n1 .eq. 1) interp(2) = interp(1)
!
    knume = '&&RFMGE1.NUME_ORDR'
    call rsutnu(modgen, ' ', 1, knume, nbordr, &
                epsi, crit, iret)
    if (iret .ne. 0) then
        call utmess('F', 'UTILITAI4_11')
    end if
    call jeveuo(knume, 'L', jordr)
!
!     --- CREATION DE LA FONCTION ---
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
    zk24(lpro) = 'FONCTION        '
    zk24(lpro+1) = interp(1)//interp(2)
    zk24(lpro+2) = 'FREQ            '
    zk24(lpro+4) = 'EE              '
    zk24(lpro+5) = nomfon
!
    call wkvect(nomfon//'.VALE', 'G V R', 2*nbordr, lvar)
    lfon = lvar+nbordr-1
!
    call getvtx(' ', 'NOM_PARA_RESU', scal=npara, nbret=n1)
    if (n1 .ne. 0) then
        zk24(lpro+3) = npara
        do iord = 1, nbordr
            call rsadpa(modgen, 'L', 1, 'FREQ', zi(jordr+iord-1), &
                        0, sjv=iad, styp=k8b)
            zr(lvar-1+iord) = zr(iad)
            call rsadpa(modgen, 'L', 1, npara, zi(jordr+iord-1), &
                        0, sjv=iad, styp=k8b)
            zr(lfon-1+iord) = zr(iad)
        end do
        goto 999
    end if
!
    call getvtx(' ', 'NOM_CHAM', scal=nomcha, nbret=n1)
    call getvis(' ', 'NUME_CMP_GENE', scal=ncmp, nbret=n1)
!
    zk24(lpro+3) = nomcha
!
    do iord = 1, nbordr
!
        call rsexch('F', modgen, nomcha, zi(jordr+iord-1), noch19, &
                    iret)
        call rsadpa(modgen, 'L', 1, 'FREQ', zi(jordr+iord-1), &
                    0, sjv=iad, styp=k8b)
        zr(lvar+iord) = zr(iad)
!
        call jeveuo(noch19//'.VALE', 'L', vr=vale)
        call jelira(noch19//'.VALE', 'TYPE', cval=k16b)
        if (k16b(1:1) .ne. 'R') then
            call utmess('F', 'UTILITAI4_17')
        end if
!
        call jeveuo(noch19//'.REFE', 'L', vk24=refe)
        mode = refe(1) (1:8)
        if (mode .eq. '        ') then
            nugene = refe(2) (1:14)
            call jeveuo(nugene//'.NUME.DEEQ', 'L', vi=deeq)
            call jeveuo(nugene//'.NUME.NEQU', 'L', vi=nequ)
            nbmode = nequ(1)
            im = 0
            do i = 1, nbmode
                istru = deeq(1+2*(i-1)+2-1)
                if (istru .lt. 0) goto 110
                im = im+1
                if (im .eq. ncmp) goto 114
110             continue
            end do
            call utmess('F', 'UTILITAI4_14')
114         continue
            im = i
        else
            im = ncmp
        end if
!
        zr(lfon+iord) = vale(im)
!
    end do
!
999 continue
!
    call jedetr(knume)
!
    call jedema()
end subroutine
