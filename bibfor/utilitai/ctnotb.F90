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
subroutine ctnotb(nbno, mesnoe, noma, nbval, nkcha, &
                  nkcmp, toucmp, nbcmp, typac, &
                  nrval, resu, nomtb, nsymb, nival, &
                  niord, label)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cnocns.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbajli.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8), parameter :: ndim = 3
    integer(kind=8) :: nbcmp, nbno, nbval
    character(len=8) :: typac, noma, resu, nomtb
    character(len=16) :: nsymb
    character(len=24) :: nkcha, nkcmp, mesnoe, nival, nrval, niord
    aster_logical :: toucmp
    character(len=24) :: label, slabel
!     ----- OPERATEUR CREA_TABLE , MOT-CLE FACTEUR RESU   --------------
!
!        BUT : REMPLISSAGE DE LA TABLE POUR UN CHAM_NO
!
!        IN     : NKCHA (K24)  : OBJET DES NOMS DE CHAMP
!                 RESU  (K8)   : NOM DU RESULTAT (SI RESULTAT,SINON ' ')
!                 NKCMP  (K24) : OBJET DES NOMS DE COMPOSANTES
!                 TOUCMP (L)   : INDIQUE SI TOUT_CMP EST RENSEIGNE
!                 NBCMP (I)    : NOMBRE DE COMPOSANTES LORSQUE
!                                NOM_CMP EST RENSEIGNE, 0 SINON
!                 TYPAC (K8)   : ACCES (ORDRE,MODE,FREQ,INST)
!                 NBVAL (I)    : NOMBRE DE VALEURS D'ACCES
!                 NOMA   (K8)  : NOM DU MAILLAGE
!                 MESNOE (K24) : OBJET DES NOMS DE NOEUD
!                 NRVAL (K16)  : OBJET DES VALEURS D'ACCES (REELS)
!                 NIVAL (K16)  : OBJET DES VALEURS D'ACCES (ENTIERS)
!                 NIORD (K16)  : NOM D'OBJET DES NUMEROS D'ORDRE
!                 NSYMB (K16)  : NOM SYMBOLIQUE DU CHAMP
!                 NBNO  (I)    : NOMBRE DE NOEUDS UTILISATEURS
!
!        IN/OUT : NOMTB (K24)  : OBJET TABLE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: jcmp, jkcha, jlno, jrval, jival, jniord, i, nbnox
    integer(kind=8) :: jcnsl, nbcmpx, n, ino, indno
    integer(kind=8) :: kcp, icmp, indcmp, ni, nr, nk, kk, nbpara
    integer(kind=8) :: j, ibid
    complex(kind=8) :: cbid
    character(len=8) :: kno
    character(len=19) :: chamns
    character(len=8), pointer :: nom_cmp(:) => null()
    character(len=16), pointer :: table_parak(:) => null()
    integer(kind=8), pointer :: table_vali(:) => null()
    character(len=16), pointer :: table_valk(:) => null()
    real(kind=8), pointer :: table_valr(:) => null()
    real(kind=8), pointer :: val_cmp(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- 0. INITIALISATIONS
    ibid = 0
    cbid = (0.d0, 0.d0)
    chamns = '&&CTNOTB.CNS       '
    call jeveuo(nkcmp, 'L', jcmp)
    call jeveuo(nkcha, 'L', jkcha)
    call jeveuo(mesnoe, 'L', jlno)
    call jeveuo(nrval, 'L', jrval)
    call jeveuo(nival, 'L', jival)
    call jeveuo(niord, 'L', jniord)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
!     TABLEAU DES VALEURS ENTIERES DE LA TABLE: ZI(JI)
!     TABLEAU DES VALEURS REELES DE LA TABLE: ZR(JR)
!     TABLEAU DES VALEURS CARACTERES DE LA TABLE: ZK16(JK)
!     POUR DES RAISONS DE PERF, CES TABLEAUX ONT ETE SORTIS DE
!     LA BOUCLE, D'OU DES DIMENSIONS EN DUR (NOMBRE SUFFISANT)
    AS_ALLOCATE(vr=table_valr, size=50)
    AS_ALLOCATE(vi=table_vali, size=50)
    AS_ALLOCATE(vk16=table_valk, size=50)

    do i = 1, nbval
!
        if (zk24(jkcha+i-1) (1:18) .ne. '&&CHAMP_INEXISTANT') then
!            -- PASSAGE CHAM_NO => CHAM_NO_S
            call cnocns(zk24(jkcha+i-1), 'V', chamns)
            call jeveuo(chamns//'.CNSV', 'L', vr=cnsv)
            call jeveuo(chamns//'.CNSL', 'L', jcnsl)
            call jeveuo(chamns//'.CNSD', 'L', vi=cnsd)
            call jeveuo(chamns//'.CNSC', 'L', vk8=cnsc)
!
!           Colonne CHAM_GD/RESULTAT
            slabel = label
            if (slabel .eq. ' ') then
                if (resu .ne. ' ') then
                    slabel = resu
                else
                    slabel = zk24(jkcha+i-1)
                end if
            end if

!             NOMBRE DE NOEUDS MAX DU CHAMP: NBNOX
            nbnox = cnsd(1)
!             NOMBRE DE COMPOSANTES MAX DU CHAMP : NBCMPX
            nbcmpx = cnsd(2)
!
!             NOMBRE DE COMPOSANTES DESIREES : N
            if (toucmp) then
                n = nbcmpx
            else
                n = nbcmp
            end if
!
!             TABLEAU DES VALEURS DES COMPOSANTES DESIREES: ZR(JVAL)
            AS_ALLOCATE(vr=val_cmp, size=n)
!
!             TABLEAU DES NOMS DE COMPOSANTES DESIREES : ZK8(JKVAL)
            AS_ALLOCATE(vk8=nom_cmp, size=n)
!
!             -- ON PARCOURT LES NOEUDS MAX,
            do ino = 1, nbnox
!
!               - SI LE NOEUD FAIT PARTIE DES NOEUDS DESIRES,
!               ON POURSUIT, SINON ON VA AU NOEUD SUIVANT:
                indno = indiis(zi(jlno), ino, 1, nbno)
                if (indno .eq. 0) goto 110
                kcp = 0
!
!              - ON PARCOURT LES COMPOSANTES:
                do icmp = 1, nbcmpx
                    if (.not. toucmp) then
!                    -SI LA COMPOSANTE FAIT PARTIE DES COMPOSANTES
!                     DESIREES, ON POURSUIT, SINON ON VA A LA
!                     COMPOSANTE SUIVANTE
                        indcmp = indik8(zk8(jcmp), cnsc(icmp), &
                                        1, nbcmp)
                        if (indcmp .eq. 0) goto 120
                    end if
!                  - SI LE CHAMP A UNE VALEUR, ON POURSUIT ET ON
!                  STOCKE LE NOM ET LA VALEUR DE COMPOSANTE :
                    if (.not. zl(jcnsl+nbcmpx*(ino-1)+icmp-1)) goto 120
                    kcp = kcp+1
                    val_cmp(kcp) = cnsv(1+nbcmpx*(ino-1)+icmp-1)
                    nom_cmp(kcp) = cnsc(icmp)
120                 continue
                end do
!
!               SOIT NI LE NOMBRE DE VALEURS ENTIERES DE LA TABLE
!               SOIT NR LE NOMBRE DE VALEURS REELES DE LA TABLE
!               SOIT NK LE NOMBRE DE VALEURS CARACTERES DE LA TABLE
!
                nr = ndim+kcp
                ni = 1
                nk = 3
                if (resu .ne. ' ') then
                    if (typac .eq. 'FREQ' .or. typac .eq. 'INST') then
                        nr = nr+1
                    else if (typac .eq. 'MODE') then
                        ni = ni+1
                    end if
                else
                    ni = 0
                    nk = 2
                end if
!
!
!               ON REMPLIT LES TABLEAUX ZI(JI),ZR(JR),ZK16(JK)
                kk = 0
                if (typac .eq. 'FREQ' .or. typac .eq. 'INST') then
                    table_valr(kk+1) = zr(jrval+i-1)
                    kk = kk+1
                end if
                do j = 1, ndim
                    table_valr(kk+1) = vale(1+3*(ino-1)+j-1)
                    kk = kk+1
                end do
                do j = 1, kcp
                    table_valr(kk+1) = val_cmp(j)
                    kk = kk+1
                end do
!
                kk = 0
                table_valk(kk+1) = slabel(1:16)
                kk = kk+1
                if (resu .ne. ' ') then
                    table_valk(kk+1) = nsymb
                    kk = kk+1
                    table_vali(1) = zi(jniord+i-1)
                    if (typac .eq. 'MODE') table_vali(1+1) = zi(jival+i-1)
                end if
                kno = int_to_char8(ino)
                table_valk(kk+1) = kno
!
!
!               TABLEAU DES NOMS DE PARAMETRES DE LA TABLE: ZK16(JPARAK)
                nbpara = nr+ni+nk
                AS_ALLOCATE(vk16=table_parak, size=nbpara)
!
!               ON REMPLIT ZK16(JPARAK)
                kk = 0
                if (resu .eq. ' ') then
                    table_parak(kk+1) = 'CHAM_GD'
                    kk = kk+1
                else
                    table_parak(kk+1) = 'RESULTAT'
                    kk = kk+1
                    table_parak(kk+1) = 'NOM_CHAM'
                    kk = kk+1
                    if (typac .ne. 'ORDRE') then
                        table_parak(kk+1) = typac
                        kk = kk+1
                    end if
                    table_parak(kk+1) = 'NUME_ORDRE'
                    kk = kk+1
                end if
                table_parak(kk+1) = 'NOEUD'
                kk = kk+1
                table_parak(kk+1) = 'COOR_X'
                kk = kk+1
                if (ndim .ge. 2) then
                    table_parak(kk+1) = 'COOR_Y'
                    kk = kk+1
                end if
                if (ndim .eq. 3) then
                    table_parak(kk+1) = 'COOR_Z'
                    kk = kk+1
                end if
                do j = 1, kcp
                    table_parak(kk+1) = nom_cmp(j)
                    kk = kk+1
                end do
!
!
!               ON AJOUTE LA LIGNE A LA TABLE
                if (resu .eq. ' ') then
                    call tbajli(nomtb, nbpara, table_parak, [ibid], table_valr, &
                                [cbid], table_valk, 0)
                else
                    call tbajli(nomtb, nbpara, table_parak, table_vali, table_valr, &
                                [cbid], table_valk, 0)
                end if
                AS_DEALLOCATE(vk16=table_parak)
!
110             continue
            end do
            AS_DEALLOCATE(vr=val_cmp)
            AS_DEALLOCATE(vk8=nom_cmp)
!
        end if
!
    end do
!
    AS_DEALLOCATE(vr=table_valr)
    AS_DEALLOCATE(vi=table_vali)
    AS_DEALLOCATE(vk16=table_valk)
!
!
    call jedema()
!
end subroutine
