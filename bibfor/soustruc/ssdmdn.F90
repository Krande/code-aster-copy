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

subroutine ssdmdn(mag)
    implicit none
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getltx.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: mag
! ----------------------------------------------------------------------
!     BUT:
!        - TRAITER LE MOT CLEF "DEFI_NOEUD"
!          DE LA COMMANDE DEFI_MAILLAGE.
!        - CREER LES OBJETS :
!            BASE VOLATILE: .NOMNOE_2
!
!     IN:
!        MAG : NOM DU MAILLAGE QUE L'ON DEFINIT.
!
    character(len=8) :: nomacr, nosma, kbid, mal, pref, nomnol, nomnog
    integer(kind=8) :: indi(4)
    character(len=24) :: valk(2)
! ----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1
    integer(kind=8) ::  ianon2, iasupm, ino, ino1
    integer(kind=8) :: inol, iocc, isma, kk, lmail, lnoeu, longt
    integer(kind=8) :: lpref, n1, n2, n3, nbnoe, nbnoet, nbnoex, lpr(1)
    integer(kind=8) :: nbnol, nbsma, nnnoe, nocc, exinno, ier
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: dime_2(:) => null()
    integer(kind=8), pointer :: conx(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: lino(:) => null()
    integer(kind=8), pointer :: noeud_conf(:) => null()
    aster_logical :: lcolle
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(mag//'.DIME', 'L', vi=dime)
    call jeveuo(mag//'.DIME_2', 'L', vi=dime_2)
    call jeveuo(mag//'.NOEUD_CONF', 'L', vi=noeud_conf)
    call jeveuo(mag//'.NOMACR', 'L', vk8=vnomacr)
    nnnoe = dime(1)
    nbsma = dime(4)
!
    call wkvect(mag//'.NOMNOE_2', 'V V K8', nnnoe, ianon2)
!
!
!     -- BOUCLE SUR LES OCCURENCES DU MOT-CLEF:
!     -----------------------------------------
    call getfac('DEFI_NOEUD', nocc)
    do iocc = 1, nocc
        call getvtx('DEFI_NOEUD', 'TOUT', iocc=iocc, scal=kbid, nbret=n1)
        if (n1 .eq. 1) then
!
!           -- CAS : TOUT: 'OUI'
!           --------------------
            lpr = 0
            call getltx('DEFI_NOEUD', 'PREFIXE', iocc, 8, 1, &
                        lpr, n2)
            lpref = lpr(1)
            call getvis('DEFI_NOEUD', 'INDEX', iocc=iocc, nbval=4, vect=indi, &
                        nbret=n3)
            lmail = indi(2)-indi(1)+1
            lnoeu = indi(4)-indi(3)+1
            lmail = max(lmail, 0)
            lnoeu = max(lnoeu, 0)
            longt = lpref+lmail+lnoeu
            if (longt .gt. 8) then
                call utmess('F', 'SOUSTRUC_57')
            end if
            if (lpref .gt. 0) then
                call getvtx('DEFI_NOEUD', 'PREFIXE', iocc=iocc, scal=pref, nbret=n2)
            end if
!
            do isma = 1, nbsma
                call jeveuo(jexnum(mag//'.SUPMAIL', isma), 'L', iasupm)
                call jenuno(jexnum(mag//'.SUPMAIL', isma), nosma)
                nomacr = vnomacr(isma)
                call jeveuo(nomacr//'.CONX', 'L', vi=conx)
                call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
                call jeexin(mal//'.NOMNOE', exinno)
                nbnoe = dime_2(4*(isma-1)+1)
                nbnol = dime_2(4*(isma-1)+2)
                nbnoet = nbnoe+nbnol
                lcolle = .false.
                call jeexin(mal//'.NOMNOE', ier)
                if (ier .ne. 0) then
                    lcolle = .true.
                end if
                do i = 1, nbnoet
                    ino = zi(iasupm-1+i)
                    if (ino .gt. nnnoe) goto 3
                    ino1 = conx(3*(i-1)+2)
                    nomnol = int_to_char8(ino1, lcolle, mal, "NOEUD")
                    if (exinno .eq. 0) then
                        nomnol = "N"//nomnol(1:7)
                    end if
                    i1 = 1
                    if (lpref .gt. 0) zk8(ianon2-1+ino) (i1:i1-1+lpref) = pref(1:lpref)
                    i1 = i1+lpref
                    if (lmail .gt. 0) zk8(ianon2-1+ino) (i1:i1-1+lmail) = nosma(indi(1):indi(2))
                    i1 = i1+lmail
                    if (lnoeu .gt. 0) zk8(ianon2-1+ino) (i1:i1-1+lnoeu) = nomnol(indi(3):indi(4))
3                   continue
                end do
            end do
        else
!
!
!       -- CAS : MAILLE, NOEUD_FIN, NOEUD_INIT :
!       ---------------------------------------
            call getvtx('DEFI_NOEUD', 'SUPER_MAILLE', iocc=iocc, scal=nosma, nbret=n1)
            call getvtx('DEFI_NOEUD', 'NOEUD_FIN', iocc=iocc, scal=nomnog, nbret=n2)
            call getvtx('DEFI_NOEUD', 'NOEUD_INIT', iocc=iocc, scal=nomnol, nbret=n3)
            if ((n1*n2*n3) .eq. 0) then
                call utmess('F', 'SOUSTRUC_58')
            end if
!
            call jenonu(jexnom(mag//'.SUPMAIL', nosma), isma)
            nomacr = vnomacr(isma)
            call jeveuo(nomacr//'.LINO', 'L', vi=lino)
            call jelira(nomacr//'.LINO', 'LONUTI', nbnoex)
            call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
            lcolle = .false.
            call jeexin(mal//'.NOMNOE', ier)
            if (ier .ne. 0) then
                lcolle = .true.
            end if
            inol = char8_to_int(nomnol, lcolle, mal, 'NOEUD')
            kk = indiis(lino, inol, 1, nbnoex)
            if (kk .eq. 0) then
                valk(1) = nomnol
                valk(2) = nosma
                call utmess('A', 'SOUSTRUC_59', nk=2, valk=valk)
                goto 1
            end if
!
            ino = dime_2(4*(isma-1)+3)+kk
            if (noeud_conf(ino) .eq. ino) then
                zk8(ianon2-1+ino) = nomnog
            else
                valk(1) = nomnol
                valk(2) = nosma
                call utmess('A', 'SOUSTRUC_60', nk=2, valk=valk)
            end if
        end if
1       continue
    end do
!
    call jedema()
end subroutine
