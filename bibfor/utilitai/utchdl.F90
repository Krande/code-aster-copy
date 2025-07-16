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

subroutine utchdl(cham19, nomma, nomail, nonoeu, nupo, &
                  nusp, ivari, nocmp1, iddl, nogranz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/numel2.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nupo, ivari, iddl, nusp
    character(len=*) :: cham19, nomma, nomail, nonoeu, nocmp1
    aster_logical, intent(in), optional :: nogranz
! ----------------------------------------------------------------------
!
! person_in_charge: jacques.pellet at edf.fr
! A_UTILI
!
! BUT: RECUPERER UN NUMERO DE DDL DANS UN CHAM_ELEM
! ----------------------------------------------------------------------
! IN  : CHAM19 : NOM DU CHAM_ELEM
! IN  : NOMMA  : NOM DU MAILLAGE
! IN  : NOMAIL : NOM DE LA MAILLE A EXTRAIRE
! IN  : NONOEU : NOM D'UN NOEUD (POUR LES CHAM_ELEM "AUX NOEUDS").
! IN  : NUPO   : NUMERO D'UN POINT (POUR LES CHAM_ELEM "GAUSS").
! IN  : NUSP   : NUMERO DU SOUS_POINT A TESTER SUR LA MAILLE NOMAIL
!                (SI NUSP=0 : IL N'Y A PAS DE SOUS-POINT)
! IN  : NOCMP1 : NOM DU DDL A EXTRAIRE SUR LE POINT CHERCHE
! IN  : IVARI  : NUMERO DE LA CMP (POUR VARI_R)
! IN  : NOGRAN : COMPORTEMENT SI LA COMPOSANTE N'EST PAS TROUVEE
!                (ARREUR FATALE SI FALSE)
! OUT : IDDL   : NUMERO DU DDL DANS LE .CELV
!   CONVENTION : IDDL=0 -> ON N'A PAS TROUVE LE DDL CHERCHE
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: ibid, gd, incmp
    integer(kind=8) :: vali(2)
    integer(kind=8) :: nec, icmp, ncmpmx, iancmp, ima
    integer(kind=8) :: ino, iaconx, nbno, ipo, nupo2, igr, iel
    integer(kind=8) :: imolo, jmolo, ispt, jlpt, nbpt, ipt, ico
    integer(kind=8) :: k, iadg, kcmp, cumu, nbspt, adiel, lgcata, ncdyn
    character(len=1) :: aof
    character(len=24) :: valk(2)
    character(len=8) :: k8b, nocmp, nomaiz, nonoez, nommaz, nomgd
    character(len=16) :: nomcmd
    character(len=19) :: noligr, chm19z, ncmp
    aster_logical :: diff, trouve, nogran, l_parallel_mesh
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: celk(:) => null()
    integer(kind=8), pointer :: long_pt_cumu(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    iddl = 0
!
    call getres(k8b, k8b, nomcmd)
    if (nomcmd .eq. 'TEST_RESU') then
        aof = 'A'
    else
        aof = 'F'
    end if
!
    chm19z = cham19
    nomaiz = nomail(1:8)
    nonoez = nonoeu(1:8)
    nommaz = nomma(1:8)
    ncmp = '&&UTCHDL.N_CMP'
    call jeveuo(chm19z//'.CELK', 'L', vk24=celk)
    noligr = celk(1) (1:19)
    trouve = .false.
    nogran = ASTER_FALSE
    if (present(nogranz)) nogran = nogranz
    l_parallel_mesh = isParallelMesh(nommaz)
!
!
!
!     1. ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
!     ---------------------------------------------------------
!     CALL CELVER(CHM19Z,'NBSPT_1','STOP',IBID)
!
!
    call jeveuo(chm19z//'.CELD', 'L', vi=celd)
    gd = celd(1)
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nomgd)
    nec = nbec(gd)
!
!
!     2. ON RECHERCHE LE NUMERO DE LA CMP CHERCHEE : ICMP
!     -------------------------------------------------------
    nocmp = nocmp1
    if (nomgd .eq. 'VARI_R') then
        nocmp = 'VARI'
        icmp = ivari
        if (icmp .le. 0) icmp = 0
    else
        call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
        call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iancmp)
        icmp = indik8(zk8(iancmp), nocmp, 1, ncmpmx)
        call wkvect(ncmp, 'V V K8', ncmpmx, incmp)
!
    end if
    if (icmp .eq. 0) then
        valk(1) = nocmp
        valk(2) = nomgd
        call utmess(aof, 'UTILITAI5_30', nk=2, valk=valk)
        iddl = 0
        goto 999
    end if
!
!
!     3. ON VERIFIE LA MAILLE : IMA
!     -----------------------------
    ima = char8_to_int(nomaiz)
    if (ima .le. 0) then
        if (.not. l_parallel_mesh) then
            valk(1) = nomaiz
            valk(2) = nommaz
            call utmess(aof, 'UTILITAI5_31', nk=2, valk=valk)
        end if
        iddl = 0
        goto 999
    end if
!
!
!     4. ON VERIFIE LE NOEUD : CALCUL DE NUPO2
!     ------------------------------------------
    if (nonoeu(1:1) .ne. ' ') then
        call dismoi('TYPE_CHAMP', chm19z, 'CHAMP', repk=k8b, arret='C', &
                    ier=ibid)
        if (k8b(1:4) .ne. 'ELNO') then
            call utmess(aof, 'UTILITAI5_32', sk=chm19z)
        end if
        ino = char8_to_int(nonoez)
        if (ino .le. 0) then
            valk(1) = nonoez
            valk(2) = nommaz
            call utmess(aof, 'ALGORITH_21', nk=2, valk=valk)
        end if
!        -- ON CHERCHE LE "IPO" CORRESPONDANT A INO:
        call jeveuo(jexnum(nommaz//'.CONNEX', ima), 'L', iaconx)
        call jelira(jexnum(nommaz//'.CONNEX', ima), 'LONMAX', nbno)
        ipo = indiis(zi(iaconx), ino, 1, nbno)
        if (ipo .le. 0) then
            valk(1) = nonoez
            valk(2) = nomaiz
            call utmess(aof, 'SOUSTRUC_59', nk=2, valk=valk)
        end if
        nupo2 = ipo
    else
        nupo2 = nupo
    end if
!
!
!     5. CALCUL DE IGR ET IEL :
!     ------------------------------------------
    call numel2(chm19z, ima, igr, iel)
    if ((igr .le. 0) .or. (iel .le. 0)) then
        valk(1) = nomaiz
        valk(2) = noligr
        call utmess(aof, 'UTILITAI5_34', nk=2, valk=valk)
    end if
    nbspt = celd(celd(4+igr)+4+4*(iel-1)+1)
    adiel = celd(celd(4+igr)+4+4*(iel-1)+4)
    ncdyn = celd(celd(4+igr)+4+4*(iel-1)+2)
!
!
!
!     6. CALCUL DE IDDL :
!     ------------------------------------------
    imolo = celd(celd(4+igr)+2)
    if (imolo .le. 0) then
        call utmess(aof, 'UTILITAI5_35', sk=nomaiz)
    end if
    call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', jmolo)
!
!
!     -- NUMERO DE SOUS_POINT :
!     ------------------------
    if (nusp .eq. 0) then
        ispt = 1
    else
        ispt = nusp
    end if
    if (ispt .gt. nbspt) then
        call utmess(aof, 'UTILITAI5_36')
        iddl = 0
        goto 999
    end if
    if ((nusp .eq. 0) .and. (nbspt .gt. 1)) then
        call utmess(aof, 'CALCULEL_1', si=nbspt)
        iddl = 0
        goto 999
    end if
!
!
!     6.1 CAS : NOMGD /= VARI_R :
!     ----------------------------
    if (nomgd .ne. 'VARI_R') then
        diff = (zi(jmolo-1+4) .gt. 10000)
        nbpt = mod(zi(jmolo-1+4), 10000)
!
!       -- SI CHAM_ELEM / ELEM (NBPT=1) ET QUE L'ON NA PAS PRECISE NUPO
!          ON PREND NUPO2=1
        if ((nupo2 .eq. 0) .and. (nbpt .eq. 1)) nupo2 = 1
!
        if (nupo2 .gt. nbpt) then
            vali(1) = nupo2
            vali(2) = nbpt
            call utmess(aof, 'UTILITAI5_37', ni=2, vali=vali)
            iddl = 0
            goto 999
        end if
        call wkvect('&&UTCHDL.LONG_PT', 'V V I', nbpt, jlpt)
        AS_ALLOCATE(vi=long_pt_cumu, size=nbpt)
!
!         -- CALCUL DU NOMBRE DE CMPS POUR CHAQUE POINT
!            ET DU CUMUL SUR LES POINTS PRECEDENTS :
        do ipt = 1, nbpt
            ico = 0
            k = 1
            if (diff) k = ipt
            iadg = jmolo-1+4+(k-1)*nec+1
            do kcmp = 1, ncmpmx
                if (exisdg(zi(iadg), kcmp)) then
                    ico = ico+1
                    zk8(incmp+ico-1) = zk8(iancmp+kcmp-1)
                    if (nocmp .eq. zk8(incmp+ico-1)) trouve = .true.
                end if
            end do
            zi(jlpt-1+ipt) = ico
        end do
        if ((.not. trouve) .and. (.not. nogran)) then
            call utmess(aof, 'UTILITAI5_38', sk=nocmp)
        end if
!
!
!
        cumu = 0
        do ipt = 1, nbpt
            long_pt_cumu(ipt) = cumu
            cumu = cumu+zi(jlpt-1+ipt)
        end do
!
!
        do ipt = 1, nbpt
            k = 1
            if (diff) k = ipt
            iadg = jmolo-1+4+(k-1)*nec+1
            ico = 0
            do kcmp = 1, ncmpmx
                if (exisdg(zi(iadg), kcmp)) then
                    ico = ico+1
!
!
                    iddl = adiel-1+nbspt*long_pt_cumu(ipt)+(ispt-1)*zi(jlpt-1+ipt)+ico
                    if ((ipt .eq. nupo2) .and. (kcmp .eq. icmp)) goto 60
!
                end if
            end do
        end do
!       -- on n'a pas trouve le point le sous-point ou la composante :
        iddl = 0
        goto 999
60      continue
!
!
!   6.2 CAS : NOMGD = VARI_R :
!   ----------------------------
    else
        lgcata = celd(celd(4+igr)+3)
        ASSERT(zi(jmolo-1+4) .le. 10000)
        nbpt = mod(zi(jmolo-1+4), 10000)
        ASSERT(nbpt .eq. lgcata)
!
        ipt = nupo2
!
        if (icmp .gt. ncdyn) then
            valk(1) = nomaiz
            vali(1) = ncdyn
            vali(2) = icmp
            call utmess(aof, 'UTILITAI7_5', sk=valk(1), ni=2, vali=vali)
            iddl = 0
            goto 999
        else
!
            if ((ispt .le. nbspt) .and. (ipt .le. nbpt)) then
                iddl = adiel-1+((ipt-1)*nbspt+ispt-1)*ncdyn+icmp
            else
                ASSERT(.false.)
            end if
        end if
    end if
!
!
999 continue
    call jedetr('&&UTCHDL.LONG_PT')
    AS_DEALLOCATE(vi=long_pt_cumu)
    call jedetr('&&UTCHDL.N_CMP')
!
    call jedema()
end subroutine
