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

subroutine ssdmgn(mag)
    implicit none
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getltx.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/indiis.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utlisi.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: mag
! ----------------------------------------------------------------------
!     BUT:
!        - TRAITER LE MOT CLEF "DEFI_GROUP_NO"
!          DE LA COMMANDE DEFI_MAILLAGE.
!        - CREER LES OBJETS :
!            BASE GLOBALE : .GROUPENO
!
!     IN:
!        MAG : NOM DU MAILLAGE QUE L'ON DEFINIT.
!
    character(len=8) :: nomacr, nomail, kbid, mal, nosma, pref
    character(len=24) :: nomgnl, nomgng
    integer(kind=8) :: indi(4)
    aster_logical :: unaun
! ----------------------------------------------------------------------
    character(len=24) :: valk(2)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i1noe, iagnl, iagno, ialino
    integer(kind=8) :: igno, ii, inol
    integer(kind=8) :: iocc, iret, isma, kk, lgnl, lmail, longt
    integer(kind=8) :: lont, lpref, n, n1, n2, n3, nbgno, lpr(1)
    integer(kind=8) :: nbgno2, nbgnot, nbid, nbno, nbnoex, nbsma, nocc
    integer(kind=8) :: nusma
    integer(kind=8), pointer :: work1(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: dime_2(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(mag//'.DIME', 'L', vi=dime)
    call jeveuo(mag//'.DIME_2', 'L', vi=dime_2)
    call jeveuo(mag//'.NOMACR', 'L', vk8=vnomacr)
    nbsma = dime(4)
    nomgng = ' '
!
!
!
!     --1 (SUR)DIMENSIONNEMENT :
!     --------------------------
    call getfac('DEFI_GROUP_NO', nocc)
!
    nbgnot = 0
    lont = 0
    do iocc = 1, nocc
        call getvis('DEFI_GROUP_NO', 'INDEX', iocc=iocc, nbval=4, vect=indi, &
                    nbret=n1)
        if (n1 .eq. 4) then
            unaun = .false.
        else
            call getvtx('DEFI_GROUP_NO', 'GROUP_NO_FIN', iocc=iocc, scal=kbid, nbret=n2)
            ASSERT(n2 .ne. 0)
            unaun = .true.
        end if
!
!
!     --1.1 CAS : INDEX, TOUT OU MAILLE :
!     -----------------------------------
        if (.not. unaun) then
            call getvtx('DEFI_GROUP_NO', 'TOUT', iocc=iocc, scal=kbid, nbret=n1)
            call getvtx('DEFI_GROUP_NO', 'SUPER_MAILLE', iocc=iocc, scal=nosma, nbret=n2)
            if (n2 .eq. 1) then
                call jeexin(jexnom(mag//'.SUPMAIL', nosma), iret)
                if (iret .eq. 0) then
                    valk(1) = nosma
                    valk(2) = mag
                    call utmess('F', 'SOUSTRUC_26', nk=2, valk=valk)
                end if
                call jenonu(jexnom(mag//'.SUPMAIL', nosma), nusma)
            end if
!
            do isma = 1, nbsma
                if ((n2 .eq. 1) .and. (nusma .ne. isma)) goto 21
                nomacr = vnomacr(isma)
                call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
                call jeexin(mal//'.GROUPENO', iret)
                if (iret .eq. 0) goto 21
                call jelira(mal//'.GROUPENO', 'NUTIOC', nbgno, kbid)
                nbgnot = nbgnot+nbgno
                do igno = 1, nbgno
                    call jelira(jexnum(mal//'.GROUPENO', igno), 'LONMAX', n3)
                    lont = lont+n3
                end do
21              continue
            end do
!
!
!     --3.2 CAS : MAILLE, GROUP_NO_FIN, GROUP_NO_INIT:
!     -----------------------------------------------
        else
            call getvtx('DEFI_GROUP_NO', 'SUPER_MAILLE', iocc=iocc, scal=nosma, nbret=n1)
            call getvtx('DEFI_GROUP_NO', 'GROUP_NO_INIT', iocc=iocc, scal=nomgnl, nbret=n)
!
            call jenonu(jexnom(mag//'.SUPMAIL', nosma), isma)
            nomacr = vnomacr(isma)
            call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
            call jelira(jexnom(mal//'.GROUPENO', nomgnl), 'LONUTI', n3)
            nbgnot = nbgnot+1
            lont = lont+n3
        end if
    end do
!
!     -- SI LONT = 0 ON S'ARRETE LA. (PLUS BAS : BUG ????)
    if ((lont .eq. 0) .or. (nbgnot .eq. 0)) goto 999
!
!
!
!
!     --2 ALLOCATION:
!     ---------------
!     --ON SURDIMENSIONNE LE NOMBRE MAX D'OBJETS DE LA COLLECTION
!       DISPERSEE .GROUPENO :
    nbgno2 = 2*nbgnot+20
    call jecrec(mag//'.GROUPENO', 'G V I', 'NOM', 'DISPERSE', 'VARIABLE', &
                nbgno2)
!
    AS_ALLOCATE(vi=work1, size=lont)
!
!
!     --3 REMPLISSAGE:
!     ----------------
    do iocc = 1, nocc
        unaun = .true.
        call getvis('DEFI_GROUP_NO', 'INDEX', iocc=iocc, nbval=4, vect=indi, &
                    nbret=n1)
        if (n1 .eq. 4) unaun = .false.
!
!
!       --3.1 CAS : INDEX, TOUT OU MAILLE :
!       -----------------------------------
        if (.not. unaun) then
            call getvtx('DEFI_GROUP_NO', 'TOUT', iocc=iocc, scal=kbid, nbret=n1)
            call getvtx('DEFI_GROUP_NO', 'SUPER_MAILLE', iocc=iocc, scal=nosma, nbret=n2)
            if (n2 .eq. 1) call jenonu(jexnom(mag//'.SUPMAIL', nosma), nusma)
            lpr = 0
            call getltx('DEFI_GROUP_NO', 'PREFIXE', iocc, 8, 1, &
                        lpr, nbid)
            lpref = lpr(1)
            call getvis('DEFI_GROUP_NO', 'INDEX', iocc=iocc, nbval=4, vect=indi, &
                        nbret=n3)
            lmail = indi(2)-indi(1)+1
            lgnl = indi(4)-indi(3)+1
            lmail = max(lmail, 0)
            lgnl = max(lgnl, 0)
            longt = lpref+lmail+lgnl
            if (longt .gt. 8) then
                call utmess('F', 'SOUSTRUC_61')
            end if
            if (lpref .gt. 0) then
                call getvtx('DEFI_GROUP_NO', 'PREFIXE', iocc=iocc, scal=pref, nbret=nbid)
            end if
!
            do isma = 1, nbsma
                if ((n2 .eq. 1) .and. (nusma .ne. isma)) goto 51
                nomacr = vnomacr(isma)
                call jenuno(jexnum(mag//'.SUPMAIL', isma), nomail)
                i1noe = dime_2(4*(isma-1)+3)
                call jeveuo(nomacr//'.LINO', 'L', ialino)
                call jelira(nomacr//'.LINO', 'LONUTI', nbnoex)
                call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
                call jeexin(mal//'.GROUPENO', iret)
                if (iret .eq. 0) then
                    nbgno = 0
                else
                    call jelira(mal//'.GROUPENO', 'NUTIOC', nbgno)
                end if
                do igno = 1, nbgno
                    call jelira(jexnum(mal//'.GROUPENO', igno), 'LONMAX', n3)
                    call jeveuo(jexnum(mal//'.GROUPENO', igno), 'L', iagnl)
                    call utlisi('INTER', zi(ialino), nbnoex, zi(iagnl), n3, &
                                work1, lont, nbno)
!
!
                    if (nbno .gt. 0) then
!
!               --3.1.1 CALCUL DE NOMGNG:
!               -------------------------
                        call jenuno(jexnum(mal//'.GROUPENO', igno), nomgnl)
                        i1 = 1
                        if (lpref .gt. 0) nomgng(i1:i1-1+lpref) = pref(1:lpref)
                        i1 = i1+lpref
                        if (lmail .gt. 0) nomgng(i1:i1-1+lmail) = nomail(indi(1):indi(2))
                        i1 = i1+lmail
                        if (lgnl .gt. 0) nomgng(i1:i1-1+lgnl) = nomgnl(indi(3):indi(4))
!
!               --3.1.2 RECOPIE DES NUMEROS DE NOEUDS:
!               --------------------------------------
                        call jecroc(jexnom(mag//'.GROUPENO', nomgng))
                        call jeecra(jexnom(mag//'.GROUPENO', nomgng), 'LONMAX', nbno)
                        call jeecra(jexnom(mag//'.GROUPENO', nomgng), 'LONUTI', nbno)
                        call jeveuo(jexnom(mag//'.GROUPENO', nomgng), 'E', iagno)
                        do ii = 1, nbno
                            inol = work1(ii)
                            kk = indiis(zi(ialino), inol, 1, nbnoex)
                            ASSERT(kk .ne. 0)
                            zi(iagno-1+ii) = i1noe+kk
                        end do
                    end if
                end do
51              continue
            end do
!
!
!       --3.2 CAS : MAILLE, GROUP_NO_FIN, GROUP_NO_INIT:
!       -----------------------------------------------
        else
            call getvtx('DEFI_GROUP_NO', 'SUPER_MAILLE', iocc=iocc, scal=nosma, nbret=n1)
            call getvtx('DEFI_GROUP_NO', 'GROUP_NO_INIT', iocc=iocc, scal=nomgnl, nbret=n)
            call getvtx('DEFI_GROUP_NO', 'GROUP_NO_FIN', iocc=iocc, scal=nomgng, nbret=n)
!
            call jenonu(jexnom(mag//'.SUPMAIL', nosma), isma)
            i1noe = dime_2(4*(isma-1)+3)
            nomacr = vnomacr(isma)
            call jeveuo(nomacr//'.LINO', 'L', ialino)
            call jelira(nomacr//'.LINO', 'LONUTI', nbnoex)
            call dismoi('NOM_MAILLA', nomacr, 'MACR_ELEM_STAT', repk=mal)
            call jelira(jexnom(mal//'.GROUPENO', nomgnl), 'LONUTI', n3)
            call jeveuo(jexnom(mal//'.GROUPENO', nomgnl), 'L', iagnl)
            call utlisi('INTER', zi(ialino), nbnoex, zi(iagnl), n3, &
                        work1, lont, nbno)
!
            if (nbno .gt. 0) then
                call jecroc(jexnom(mag//'.GROUPENO', nomgng))
                call jeecra(jexnom(mag//'.GROUPENO', nomgng), 'LONMAX', nbno)
                call jeecra(jexnom(mag//'.GROUPENO', nomgng), 'LONUTI', nbno)
                call jeveuo(jexnom(mag//'.GROUPENO', nomgng), 'E', iagno)
                do ii = 1, nbno
                    inol = work1(ii)
                    kk = indiis(zi(ialino), inol, 1, nbnoex)
                    ASSERT(kk .ne. 0)
                    zi(iagno-1+ii) = i1noe+kk
                end do
            else
                call utmess('A', 'SOUSTRUC_62', sk=nomgng)
            end if
        end if
    end do
!
!
    AS_DEALLOCATE(vi=work1)
!
999 continue
    call jedema()
end subroutine
