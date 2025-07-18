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

subroutine op0174()
    implicit none
!
!     COMMANDE:  RECU_TABLE
!
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getltx.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ltnotb.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexpa.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsorac.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: iret, lonord(1), iord, ipara, i, n, numord, ibid
    integer(kind=8) :: nbpara, inom, ityp, ilong, itabi, itabr, itabc
    integer(kind=8) :: pi, pc, pr, nbac, nbpa, inom2
    real(kind=8) :: rbid
    complex(kind=8) :: cbid
    character(len=8) :: table, concpt, typ
    character(len=16) :: nomcmd, typcon, nomsym, nom
    character(len=24) :: nomtab, kbid
!
!     ------------------------------------------------------------------
!
!
    call jemarq()
    call infmaj()
!
    call getres(table, typcon, nomcmd)
    call getvid(' ', 'CO', scal=concpt, nbret=iret)
!
!
! --------------------------
!   EXTRACTION D'UNE TABLE
! --------------------------
!
    call getvtx(' ', 'NOM_TABLE', scal=nomsym, nbret=iret)
    if (iret .ne. 0) then
        call ltnotb(concpt, nomsym, nomtab)
        call copisd('TABLE', 'G', nomtab, table)
        call titre()
        goto 9999
    end if
!
!
!
! ----------------------------
!   EXTRACTION DE PARAMETRES
! ----------------------------
!

!
! -- LECTURE DES NUMEROS D'ORDRE
!
    call rsorac(concpt, 'LONUTI', 0, rbid, kbid, &
                cbid, rbid, kbid, lonord, 1, &
                ibid)
    call wkvect('&&OP0174.NUME_ORDRE', 'V V I', lonord(1), iord)
    call rsorac(concpt, 'TOUT_ORDRE', 0, rbid, kbid, &
                cbid, rbid, kbid, zi(iord), lonord(1), &
                ibid)

!    TOUS LES PARAMETRES : MOT-CLE TOUT_PARA
!    ---------------------------------------
    call getvtx(' ', 'TOUT_PARA', scal=kbid, nbret=iret)
    if (iret .ne. 0) then
        ASSERT(kbid(1:3) .eq. 'OUI')

! -- -- nombre de parametres présents
        call rsnopa(concpt, 2, '&&OP0174.NOM_PAR2', nbac, nbpa)
        nbpara = nbac+nbpa
        if (nbpara .eq. 0) call utmess('F', 'UTILITAI3_3')

! -- -- liste des parametres
        call jeexin('&&OP0174.NOM_PAR2', iret)
        ASSERT(iret .gt. 0)
        call jeveuo('&&OP0174.NOM_PAR2', 'E', jadr=inom2)

! -- -- ajout NUME_ORDRE en début de liste
        call wkvect('&&OP0174.NOM_PARA', 'V V K16', nbpara+1, inom)
        zk16(inom) = 'NUME_ORDRE'

! -- -- sélection des types i, r et c uniquement
        call wkvect('&&OP0174.TYPE_PARA', 'V V K8 ', nbpara+1, ityp)
        zk8(ityp) = 'I'
        numord = zi(iord)
        nbpa = 0
        do i = 1, nbpara
            nom = zk16(inom2-1+i)
            call rsadpa(concpt, 'L', 1, nom, numord, &
                        1, sjv=ipara, styp=typ, istop=0)
            if (typ(1:1) .eq. 'R' .or. typ(1:1) .eq. 'I' &
                .or. typ(1:1) .eq. 'C') then
                nbpa = nbpa+1
                zk16(inom+nbpa) = zk16(inom2-1+i)
                zk8(ityp+nbpa) = typ(1:1)
            end if
        end do
        nbpara = nbpa
        call jedetr('&&OP0174.NOM_PAR2')
    else

!    NOMBRE DE PARAMETRES : MOT-CLE NOM_PARA
!    ---------------------------------------
        call getvtx(' ', 'NOM_PARA', nbval=0, nbret=nbpara)
        ASSERT(nbpara .ne. 0)
        nbpara = -nbpara

! -- -- noms des paramètres à extraire
        call wkvect('&&OP0174.NOM_PARA', 'V V K16', nbpara+1, inom)
        call wkvect('&&OP0174.LONG_NOM_PARA', 'V V I', nbpara, ilong)
!
        call getltx(' ', 'NOM_PARA', 1, 16, nbpara, &
                    zi(ilong), ibid)
        do i = 1, nbpara
            if (zi(ilong-1+i) .gt. 16) then
                call utmess('F', 'UTILITAI3_4')
            end if
        end do
!
        zk16(inom) = 'NUME_ORDRE'
        call getvtx(' ', 'NOM_PARA', nbval=nbpara, vect=zk16(inom+1), nbret=ibid)
!
        do i = 1, nbpara
            nom = zk16(inom+i)
            call rsexpa(concpt, 2, nom, iret)
            if (iret .eq. 0) then
                call utmess('F', 'UTILITAI3_5', sk=nom)
            end if
        end do

! -- -- types des paramètres à extraire
        call wkvect('&&OP0174.TYPE_PARA', 'V V K8 ', nbpara+1, ityp)
        zk8(ityp) = 'I'
        numord = zi(iord)
        do i = 1, nbpara
            nom = zk16(inom+i)
            call rsadpa(concpt, 'L', 1, nom, numord, &
                        1, sjv=ipara, styp=typ, istop=0)
            if (typ(1:1) .ne. 'R' .and. typ(1:1) .ne. 'I' .and. typ(1:1) .ne. 'C') then
                call utmess('F', 'UTILITAI3_6')
            end if
            zk8(ityp+i) = typ(1:1)
        end do
    end if
!
!
! -- INITIALISATION DE LA TABLE
!
    call tbcrsd(table, 'G')
    call titre()
    call tbajpa(table, 1+nbpara, zk16(inom), zk8(ityp))
!
!
! -- EXTRACTION DES PARAMETRES
!
    call wkvect('&&OP0174.PARA_R', 'V V R', nbpara+1, itabr)
    call wkvect('&&OP0174.PARA_I', 'V V I', nbpara+1, itabi)
    call wkvect('&&OP0174.PARA_C', 'V V C', nbpara+1, itabc)
!
    do n = 1, lonord(1)
        numord = zi(iord-1+n)
!
        pi = 0
        pr = 0
        pc = 0
!
        zi(itabi+pi) = numord
        pi = pi+1
!
        do i = 1, nbpara
            nom = zk16(inom+i)
            call rsadpa(concpt, 'L', 1, nom, numord, &
                        1, sjv=ipara, styp=typ, istop=0)
            if (typ(1:1) .eq. 'R') then
                zr(itabr+pr) = zr(ipara)
                pr = pr+1
            else if (typ(1:1) .eq. 'I') then
                zi(itabi+pi) = zi(ipara)
                pi = pi+1
            else if (typ(1:1) .eq. 'C') then
                zc(itabc+pc) = zc(ipara)
                pc = pc+1
            else
                ASSERT(.false.)
            end if
        end do
!
        call tbajli(table, 1+nbpara, zk16(inom), zi(itabi), zr(itabr), &
                    zc(itabc), kbid, 0)
    end do
!
!
!
9999 continue
    call jedema()
end subroutine
