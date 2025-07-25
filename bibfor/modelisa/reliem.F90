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

subroutine reliem(mo, ma, typem, motfaz, iocc, &
                  nbmocl, limocl, tymocl, litroz, nbtrou, l_keep_propz, l_allz)
    implicit none
#include "asterf_types.h"
#include "asterfort/addPhantomNodesFromCells.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
#include "jeveux.h"
!
    integer(kind=8) :: iocc, nbmocl, nbtrou
    character(len=8) :: ma, modele
    character(len=*) :: limocl(nbmocl), tymocl(nbmocl), mo
    character(len=*) :: litroz, typem, motfaz
    aster_logical, optional, intent(in) :: l_keep_propz
    aster_logical, optional, intent(in) :: l_allz
! ----------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
!
!     CE MODULE PERMET DE CREER UN OBJET JEVEUX CONTENANT UNE LISTE
!     DE NOMS OU NUMEROS DE MAILLES OU DE NOEUDS CORRESPONDANT AUX
!     MOTS-CLES TRANSMIS EN ARGUMENTS.
!
! IN  : MO     : NOM DU MODELE (FACULTATIF SINON : ' ')
!           SI LE NOM DU MODELE EST DONNE, ON VERIFIERA QUE LES MAILLES
!           (OU LES NOEUDS) RECUPERES FONT PARTIE DU MODELE.
!           S'ILS NE FONT PAS PARTIE DU MODELE => ALARME
! IN  : MA     : NOM DU MAILLAGE
! IN  : TYPEM  : PRECISE LE TYPE DE LISTE QUE L'ON VEUT RECUPERER
!              : 'NU_MAILLE'  : NUMEROS DE MAILLES
!              : 'NO_MAILLE'  : NOMS    DE MAILLES
!              : 'NU_NOEUD'   : NUMEROS DE NOEUDS
!              : 'NO_NOEUD'   : NOMS    DE NOEUDS
! IN  : MOTFAZ : NOM DU MOT CLE FACTEUR (OU ' ')
! IN  : IOCC   : NUMERO DE L'OCCURENCE DU MOT CLE FACTEUR
! IN  : NBMOCL : NOMBRE DE MOTS CLES A SCRUTER
!                (DIMENSION DE LIMOCL)
! IN  : LIMOCL : LISTE DES MOTS CLE A SCRUTER
! IN  : TYMOCL : LISTE DES TYPES DE MOTS CLE A SCRUTER :
!                / 'GROUP_MA'
!                / 'GROUP_NO'
!                / 'MAILLE'
!                / 'NOEUD'
!                / 'TOUT'   % TOUT:'OUI'
! IN/JXOUT : LITROZ : NOM DE L'OBJET JEVEUX QUI CONTIENDRA LA LISTE DES
!                     ENTITES (MAILLE OU NOEUD) TROUVEES
! OUT : NBTROU : NOMBRE D'ENTITES TROUVEES
! IN, OPTIONAL : L_KEEP_PROP : PRIS EN COMPTE UNIQUEMENT UN MAILLAGE PARALLELE
!    (CELA NE CHANGE RIEN DANS LES AUTRES CAS)
!    POUR UN PARALLEL_MESH, SI TRUE ON NE GARDE QUE LES MAILLES/NOEUDS DONT LE SOUS-DOMAINE
!    EST PROPRIETAIRE SI FALSE ON GARDE TOUT
!    SI L'ARGUMENT N'EST PAS PRESENT ON GARDE TOUT (=FALSE).
! IN, OPTIONAL     : L_ALL : TRUE  : forcer comme TOUT='OUI'
!                            FALSE : par défaut, cas normal
! ----------------------------------------------------------------------
    character(len=24) :: litrou
    integer(kind=8) :: jno, jma, kno, kma, iacnex, iem, nem, numno, nno, nma, nbenc
    integer(kind=8) :: ibid, ient, rank
    integer(kind=8) ::  itrma, ima, ino, nbma, nbno, nbnoma, imo, ier
    integer(kind=8) :: lma, lno, itbma, itbno, inoem, ntou, k, ifm, niv
    character(len=8) :: type2, oui, noent, nomgd
    character(len=16) :: motfac, motcle, typmcl, phenom
    character(len=19) :: ligrel
    character(len=24) :: karg
    mpi_int :: mrank
    aster_logical :: l_parallel_mesh, l_group_ma, l_keep_prop, l_all
    aster_logical :: lcolle, lcolle2
    integer(kind=8), pointer :: maille(:) => null()
    integer(kind=8), pointer :: prnm(:) => null()
    integer(kind=8), pointer :: v_maex(:) => null()
    integer(kind=8), pointer :: v_noex(:) => null()
    integer(kind=4), pointer :: indic_noeud(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    litrou = litroz
    motfac = motfaz
    modele = mo
    l_parallel_mesh = isParallelMesh(ma)
    if (present(l_keep_propz)) then
        l_keep_prop = l_keep_propz
    else
        l_keep_prop = ASTER_FALSE
    end if
    if (present(l_allz)) then
        l_all = l_allz
    else
        l_all = ASTER_FALSE
    end if
    call infniv(ifm, niv)
    lcolle = .false.
    call jeexin(ma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    lcolle2 = .false.
    call jeexin(ma//'.NOMMAI', ier)
    if (ier .ne. 0) then
        lcolle2 = .true.
    end if
!
    call asmpi_info(rank=mrank)
    rank = to_aster_int(mrank)
!
!     --- VERIFICATIONS PRELIMINAIRES ---
!
    if (typem .ne. 'NO_MAILLE' .and. typem .ne. 'NO_NOEUD' .and. typem .ne. 'NU_MAILLE' &
        .and. typem .ne. 'NU_NOEUD') then
        ASSERT(.false.)
    end if
!
    type2 = typem(4:)
    l_group_ma = ASTER_FALSE
    do imo = 1, nbmocl
        motcle = limocl(imo)
        typmcl = tymocl(imo)
        if (typmcl .eq. 'NOEUD' .or. typmcl .eq. 'GROUP_NO') then
            if (type2 .eq. 'MAILLE') then
                ASSERT(.false.)
            end if
        else if (typmcl .ne. 'MAILLE' .and. typmcl .ne. 'GROUP_MA' .and. typmcl .ne. 'TOUT') then
            ASSERT(.false.)
        end if
        if (typmcl == 'GROUP_MA') then
            l_group_ma = ASTER_TRUE
        end if
    end do
!
!   EN CAS D'EXISTENCE DE L'OBJET, ON LE DETRUIT
    call jedetr(litrou)
!
!     --- CREATION DES TABLEAUX DE TRAVAIL ---
!
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    if (nbma .gt. 0) then
        call wkvect('&&RELIEM.INDIC_MAILLE', 'V V S', max(nbma, 1), itrma)
        if (modele .ne. ' ') call jeveuo(modele//'.MAILLE', 'L', vi=maille)
    end if
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
    AS_ALLOCATE(vi4=indic_noeud, size=nbno)
!
    do k = 1, nbma
        zi4(itrma-1+k) = 0
    end do
    do k = 1, nbno
        indic_noeud(k) = 0
    end do
!
!
!
!     --- CONSTITUTION DES LISTES DES MAILLES ET DES NOEUDS
!         PAR MARQUAGE DANS LES TABLEAUX DE TRAVAIL         ---
!
!
    do imo = 1, nbmocl
        motcle = limocl(imo)
        typmcl = tymocl(imo)
!
!        -- CAS TOUT:'OUI'
!        -----------------
        if (typmcl .eq. 'TOUT') then
            if (l_all) then
                ntou = 1
            else
                call getvtx(motfac, motcle, iocc=iocc, scal=oui, nbret=ntou)
            end if
            if (ntou .gt. 0) then
                if (type2 .eq. 'MAILLE') then
                    do k = 1, nbma
                        if (modele .ne. ' ') then
                            if (maille(k) .ne. 0) then
                                zi4(itrma-1+k) = 1
                            end if
                        else
                            zi4(itrma-1+k) = 1
                        end if
                    end do
                end if
                if (type2 .eq. 'NOEUD') then
                    do k = 1, nbno
                        indic_noeud(k) = 1
                    end do
                end if
            end if
            goto 90
        end if
!
!
        call getvem(ma, typmcl, motfac, motcle, iocc, &
                    0, karg, nem)
        nem = -nem
        if (nem .eq. 0) goto 90
        if (typmcl(1:6) .ne. 'GROUP_') then
            call wkvect('&&RELIEM.NOM_EM', 'V V K8', nem, inoem)
            call getvem(ma, typmcl, motfac, motcle, iocc, &
                        nem, zk8(inoem), nem)
        else
            call wkvect('&&RELIEM.NOM_EM', 'V V K24', nem, inoem)
            call getvem(ma, typmcl, motfac, motcle, iocc, &
                        nem, zk24(inoem), nem)
        end if
!
        do iem = 1, nem
            if (typmcl(1:6) .ne. 'GROUP_') then
                karg = zk8(inoem-1+iem)
            else
                karg = zk24(inoem-1+iem)
            end if
!
            if (typmcl .eq. 'MAILLE') then
                if (l_parallel_mesh) then
                    call utmess('F', 'MODELISA7_86')
                end if
                ima = char8_to_int(karg, lcolle2, ma, "MAILLE")
                zi4(itrma-1+ima) = 1
!
            else if (typmcl .eq. 'GROUP_MA') then
                call jelira(jexnom(ma//'.GROUPEMA', karg), 'LONUTI', nma)
                call jeveuo(jexnom(ma//'.GROUPEMA', karg), 'L', kma)
!
!           -- UNE VERIFICATION PENDANT LE CHANTIER "GROUPES VIDES" :
                call jelira(jexnom(ma//'.GROUPEMA', karg), 'LONMAX', ibid)
                if (ibid .eq. 1) then
                    ASSERT(nma .le. 1)
                else
                    ASSERT(nma .eq. ibid)
                end if
!
                do jma = 1, nma
                    ima = zi(kma-1+jma)
                    zi4(itrma-1+ima) = 1
                end do
!
            else if (typmcl .eq. 'NOEUD') then
                if (l_parallel_mesh) then
                    call utmess('F', 'MODELISA7_86')
                end if
                ino = char8_to_int(karg, lcolle, ma, "NOEUD")
                indic_noeud(ino) = 1
!
            else if (typmcl .eq. 'GROUP_NO') then
                call jelira(jexnom(ma//'.GROUPENO', karg), 'LONUTI', nno)
                call jeveuo(jexnom(ma//'.GROUPENO', karg), 'L', kno)
!
!           -- UNE VERIFICATION PENDANT LE CHANTIER "GROUPES VIDES" :
                call jelira(jexnom(ma//'.GROUPENO', karg), 'LONMAX', ibid)
                if (ibid .eq. 1) then
                    ASSERT(nno .le. 1)
                else
                    ASSERT(nno .eq. ibid)
                end if
!
                do jno = 1, nno
                    ino = zi(kno-1+jno)
                    indic_noeud(ino) = 1
                end do
            end if
        end do
        call jedetr('&&RELIEM.NOM_EM')
90      continue
    end do
!
!     --- AJOUT DES NOEUDS DE LA LISTE DES MAILLES A CELLE DES NOEUDS
!
    if (type2 .eq. 'NOEUD') then
        do ima = 1, nbma
            if (zi4(itrma-1+ima) .ne. 0) then
                call jeveuo(jexnum(ma//'.CONNEX', ima), 'L', iacnex)
                call jelira(jexnum(ma//'.CONNEX', ima), 'LONMAX', nbnoma)
                do ino = 1, nbnoma
                    numno = zi(iacnex-1+ino)
                    indic_noeud(numno) = 1
                end do
            end if
        end do
!
!   --- Pour un maillage ParallelMesh, il faut aussi savoir si des noeuds fantome
!       font parti aussi de mailles non-présentes localement
!
        if (l_group_ma .and. l_parallel_mesh) then
            call addPhantomNodesFromCells(ma, indic_noeud)
        end if
    end if
!
!
!     --- CREATION DE L'OBJET JEVEUX LITROU ---
    if (type2 .eq. 'MAILLE') then
!
!        --- COMPTAGE DES MAILLES ---
!
        if (l_parallel_mesh) then
            call jeveuo(ma//'.MAEX', 'L', vi=v_maex)
        end if

        nbtrou = 0
        do ima = 1, nbma
            if (zi4(itrma-1+ima) .ne. 0) then
                if (l_parallel_mesh .and. l_keep_prop) then
                    if (v_maex(ima) == rank) then
                        nbtrou = nbtrou+1
                    else
                        zi4(itrma-1+ima) = 0
                    end if
                else
                    nbtrou = nbtrou+1
                end if
            end if
        end do
        if (nbtrou .eq. 0) goto 200
!
!
!
        if (typem(1:2) .eq. 'NU') then
            call wkvect(litrou, 'V V I', nbtrou, itbma)
!
!
!           --- RANGEMENT DES NUMEROS DE MAILLES ---
            lma = 0
            do ima = 1, nbma
                if (zi4(itrma-1+ima) .ne. 0) then
                    lma = lma+1
                    zi(itbma-1+lma) = ima
                end if
            end do
!
        else
            call wkvect(litrou, 'V V K8', nbtrou, itbma)
!
!
!           --- RANGEMENT DES NOMS DE MAILLES ---
            lma = 0
            do ima = 1, nbma
                if (zi4(itrma-1+ima) .ne. 0) then
                    lma = lma+1
                    zk8(itbma-1+lma) = int_to_char8(ima, lcolle2, ma, "MAILLE")
                end if
            end do
        end if
!
!
!
!       -- ON VERIFIE QUE LES MAILLES FONT PARTIE DU MODELE :
!       ----------------------------------------------------
        if (modele .ne. ' ') then
            ier = 0
            do ima = 1, nbma
                if (zi4(itrma-1+ima) .ne. 0) then
                    if (maille(ima) .eq. 0) then
                        ier = ier+1
                        noent = int_to_char8(ima, lcolle2, ma, "MAILLE")
                        write (ifm, *) ' MAILLE : ', noent
                    end if
                end if
            end do
            if (ier .ne. 0) then
                call utmess('F', 'MODELISA6_96', sk=motfac, si=ier)
            end if
        end if
!
!
!
    else
!
!        --- COMPTAGE DES NOEUDS ---
!
        if (l_parallel_mesh) then
            call jeveuo(ma//'.NOEX', 'L', vi=v_noex)
        end if
!
        nbtrou = 0
        do ino = 1, nbno
            if (indic_noeud(ino) .ne. 0) then
                if (l_parallel_mesh .and. l_keep_prop) then
                    if (v_noex(ino) == rank) then
                        nbtrou = nbtrou+1
                    else
                        indic_noeud(ino) = 0
                    end if
                else
                    nbtrou = nbtrou+1
                end if
            end if
        end do
        if (nbtrou .eq. 0) goto 200
!
!
!
        if (typem(1:2) .eq. 'NU') then
            call wkvect(litrou, 'V V I', nbtrou, itbno)
!
!
!           --- RANGEMENT DES NUMEROS DE NOEUDS ---
            lno = 0
            do ino = 1, nbno
                if (indic_noeud(ino) .ne. 0) then
                    lno = lno+1
                    zi(itbno-1+lno) = ino
                end if
            end do
!
        else
            call wkvect(litrou, 'V V K8', nbtrou, itbno)
!
!
!           --- RANGEMENT DES NOMS DE NOEUDS ---
            lno = 0
            do ino = 1, nbno
                if (indic_noeud(ino) .ne. 0) then
                    lno = lno+1
                    zk8(itbno-1+lno) = int_to_char8(ino, lcolle, ma, "NOEUD")
                end if
            end do
        end if
!
!       -- ON VERIFIE QUE LES NOEUDS FONT PARTIE DU MODELE :
!       ----------------------------------------------------
        if (modele .ne. ' ') then
            call dismoi('PHENOMENE', modele, 'MODELE', repk=phenom)
            call dismoi('NOM_GD', phenom, 'PHENOMENE', repk=nomgd)
            call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nbenc)
            call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrel)
            call jeveuo(ligrel//'.PRNM', 'L', vi=prnm)
            ier = 0
            do ino = 1, nbno
                if (indic_noeud(ino) .ne. 0) then
                    do ient = 1, nbenc
                        if (prnm(nbenc*(ino-1)+ient) .ne. 0) goto 191
                    end do
!             LE NOEUD NE PORTE AUCUNE COMPOSANTE DE LA GRANDEUR
!             ASSOCIEE AU PHENOMENE
                    ier = ier+1
                    noent = int_to_char8(ino, lcolle, ma, "NOEUD")
                    write (ifm, *) ' NOEUD : ', noent
                end if
191             continue
            end do
            if (ier .ne. 0) then
                call utmess('F', 'MODELISA6_13', sk=motfac, si=ier)
            end if
        end if
!
    end if
!
!
!     --- DESTRUCTION DES TABLEAUX DE TRAVAIL ---
200 continue
    call jedetr('&&RELIEM.INDIC_MAILLE')
    AS_DEALLOCATE(vi4=indic_noeud)
!
    call jedema()
end subroutine
