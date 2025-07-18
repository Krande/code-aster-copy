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

subroutine w039c1(carte, ifi, form, ligrel, titre)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/exisd.h"
#include "asterfort/imprsd.h"
#include "asterfort/irceme.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/w039c2.h"
#include "asterfort/w039c4.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=19) :: ligrel
    character(len=*) :: carte, titre, form
    integer(kind=8) :: ifi
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
!     BUT:
!       IMPRIMER UNE "CARTE" D'1 CONCEPT CHAM_MATER, CARA_ELE, ...
! ----------------------------------------------------------------------
!
!
!
!
    integer(kind=8) :: iret, ima, nbma, izone, nuzone
    integer(kind=8) ::  jcesd, jcesl, iad, dec1, dec2, ifm, ifr, nncp, iexi, nbCmpDyna
    integer(kind=8) :: jdesc, jvale, ngedit, nugd, ncmpmx, kgedit, kzone, kcmp
    character(len=19) :: cart1, cel2, ces2
    character(len=64) :: nommed
    character(len=8) :: kbid, ma, tsca, nomgd, modele, typech, sdcarm
    character(len=16) :: field_type
    real(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: ptma(:) => null()
    integer(kind=8), pointer :: zones(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()

!     -- si la carte n'existe pas, il n'y a rien a faire :
!     -----------------------------------------------------
    call exisd('CARTE', carte, iexi)
    if (iexi .eq. 0) goto 999

!
    ifm = iunifi('MESSAGE')
    ifr = iunifi('RESULTAT')
    cart1 = carte

!     -- POUR QUE LE CHAM_ELEM QUE L'ON VA IMPRIMER AIT UN NOM "PROCHE"
!        DE CELUI DE LA VRAIE CARTE
    cel2 = cart1
    cel2(9:9) = '_'

!     -- QUELQUES INFOS SUR LA CARTE :
    call jeveuo(cart1//'.DESC', 'L', jdesc)
    call jeveuo(cart1//'.VALE', 'L', jvale)
    ngedit = zi(jdesc-1+3)
    if (ngedit .eq. 0) goto 999

    nugd = zi(jdesc-1+1)
    call jenuno(jexnum('&CATA.GD.NOMGD', nugd), nomgd)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', ncmpmx)
    write (ifm, *) ' '
    write (ifr, *) ' '
    write (ifm, '(A,A)') 'IMPRESSION D''UN CHAMP DE CONCEPT : ', titre
    write (ifr, '(A,A)') 'IMPRESSION D''UN CHAMP DE CONCEPT : ', titre
    write (ifm, '(A,A)') 'NOM DU CHAMP : ', cel2
    write (ifr, '(A,A)') 'NOM DU CHAMP : ', cel2

!     -- SI LA CARTE A DES VALEURS REELLES ET QUE LE FORMAT EST 'MED'
!        ON L'IMPRIME AVEC SES VALEURS REELLES. C'EST PLUS JOLI !
    if (form .eq. 'MED' .and. tsca .eq. 'R') then
        call w039c4(carte, ifi, form)
        goto 999
    end if

    write (ifm, '(A)') 'CORRESPONDANCE VALEUR <-> CONTENU :'
    write (ifr, '(A)') 'CORRESPONDANCE VALEUR <-> CONTENU :'

!     -- PARFOIS LA CARTE CONTIENT DES ZONES AYANT LES MEMES VALEURS :
!     ----------------------------------------------------------------
    AS_ALLOCATE(vi=zones, size=ngedit)
    nuzone = 0
    do kgedit = 1, ngedit
        izone = kgedit
!       -- ON REGARDE SI LES VALEURS DE IZONE N'ONT PAS DEJA ETE VUES
!          POUR KZONE < IZONE :
        do kzone = 1, izone-1
            do kcmp = 1, ncmpmx
                dec1 = ncmpmx*(kzone-1)+kcmp
                dec2 = ncmpmx*(izone-1)+kcmp
                if (tsca .eq. 'K8') then
                    if (zk8(jvale-1+dec1) .ne. zk8(jvale-1+dec2)) goto 20
                else if (tsca .eq. 'K16') then
                    if (zk16(jvale-1+dec1) .ne. zk16(jvale-1+dec2)) goto 20
                else if (tsca .eq. 'K24') then
                    if (zk24(jvale-1+dec1) .ne. zk24(jvale-1+dec2)) goto 20
                else if (tsca .eq. 'I') then
                    if (zi(jvale-1+dec1) .ne. zi(jvale-1+dec2)) goto 20
                else if (tsca .eq. 'R') then
                    if (zr(jvale-1+dec1) .ne. zr(jvale-1+dec2)) goto 20
                else if (tsca .eq. 'C') then
                    if (zc(jvale-1+dec1) .ne. zc(jvale-1+dec2)) goto 20
                else
                    ASSERT(.false.)
                end if
            end do
!         -- IZONE == KZONE :
            zones(izone) = zones(kzone)
            goto 30
!
20          continue
        end do
        nuzone = nuzone+1
        zones(izone) = nuzone
        call w039c2(izone, nuzone, jvale, jdesc, nomgd, ifm, &
                    ifr)
30      continue
    end do

!     -- ON TRANSFORME LA CARTE EN UN CHAM_ELEM_S DE NEUT_R :
!     ------------------------------------------------------
    call jelira(cart1//'.DESC', 'DOCU', cval=kbid)
    ASSERT(kbid .eq. 'CART')
    call dismoi('NOM_MAILLA', cart1, 'CARTE', repk=ma)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
!
    call etenca(cart1, ligrel, iret)
    ASSERT(iret .eq. 0)
    call jeveuo(cart1//'.PTMA', 'L', vi=ptma)
!
    ces2 = '&&W039C1.CES2'
    call cescre('V', ces2, 'ELEM', ma, 'NEUT_R', &
                1, 'X1', [-1], [-1], [-1])
    call jeveuo(ces2//'.CESD', 'L', jcesd)
    call jeveuo(ces2//'.CESV', 'E', vr=cesv)
    call jeveuo(ces2//'.CESL', 'E', jcesl)
    do ima = 1, nbma
        izone = ptma(ima)
        if (izone .gt. 0) then
            nuzone = zones(izone)
            ASSERT(nuzone .gt. 0)
            call cesexi('C', jcesd, jcesl, ima, 1, &
                        1, 1, iad)
            ASSERT(iad .le. 0)
            zl(jcesl-1-iad) = .true.
            cesv(1-1-iad) = dble(nuzone)
        end if
    end do

!     -- TRANSFORMATION DE CES2 EN CEL2 (CHAM_ELEM/ELEM) :
!     ----------------------------------------------------
    call cescel(ces2, ligrel, 'TOU_INI_ELEM', 'PNEU1_R', 'OUI', &
                nncp, 'V', cel2, 'F', iret)
    ASSERT(iret .eq. 0)
    call detrsd('CHAM_ELEM_S', ces2)

!     -- IMPRESSION DE CEL2 :
!     -----------------------

    if (form .eq. 'MED') then
!     -------------------------
        nommed = cel2
        typech = 'ELEM'
        modele = ' '
        sdcarm = ' '
        field_type = 'Unknown'
        call irceme(ifi, nommed, cel2, typech, modele, &
                    0, ' ', ' ', ' ', 0, &
                    0.d0, 0, 0, [0], sdcarm, sdcarm, &
                    field_type, nbCmpDyna, .false._1, iret)
        ASSERT(iret .eq. 0)

    else if (form .eq. 'RESULTAT') then
!     ---------------------------
        call imprsd('CHAMP', cel2, ifi, titre)

    else
        ASSERT(.false.)
    end if
    call detrsd('CHAM_ELEM', cel2)
    AS_DEALLOCATE(vi=zones)

!
999 continue
    call jedema()
end subroutine
