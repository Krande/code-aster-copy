! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine pacoap(lisi1z, lisi2z, lonlis, centre, theta, &
                  t, nomaz, liso1z, liso2z)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8gaem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/matrot.h"
#include "asterfort/padist.h"
#include "asterfort/parotr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer :: lonlis
    character(len=*) :: lisi1z, lisi2z, nomaz, liso1z, liso2z
    real(kind=8) :: centre(3), theta(3), t(3)
!     BUT: TRIER 2 LISTES DE NOEUDS LISI1Z ET LISI2Z DE MANIERE A
!     METTRE EN VIS A VIS LES NOEUDS DES 2 LISTES VIA
!     LA TRANSFORMATION  OM.THETA+T.
!     LES LISTES TRIEES OBTENUES A PARTIR DE LISI1Z ET LISI2Z
!     SONT RESPECTIVEMENT LISO1Z ET LISO2Z, LA CORRESPONDANCE
!     ENTRE LES NOEUDS DES 2 LISTES EST ASSUREE DE LA MANIERE
!     SUIVANTE :
!          POUR I =1, LONLIS
!          LISO1Z(I) EST EN VIS-AVIS AVEC LISO2Z(I)
!
!     LES LISTES LISI1Z, LISI2Z, LISO1Z ET LISO2Z CONTIENNENT
!     LES NOMS DES NOEUDS (CE SONT DES LISTES DE K8).
!
!---------------------------------------------------------------------
! ARGUMENTS D'ENTREE:
! IN   LISI1Z     K24 : NOM DE LA 1ERE LISTE
! IN   LISI2Z     K24 : NOM DE LA 2EME LISTE
! IN   LONLIS     I   : LONGUEUR COMMUNE DE CES 2 LISTES
! IN   CENTRE(3)  R   : COORDONNEES DU CENTRE DE ROTATION
! IN   THETA(3)   R   : ANGLES DE ROTATION
! IN   T(3)       R   : COORDONNEES DE LA TRANSLATION
! IN   NOMAZ      K8  : NOM DU MAILLAGE
! OUT  LISO1Z     K24 : NOM DE LA 1ERE LISTE TRIEE
! OUT  LISO2Z     K24 : NOM DE LA 2EME LISTE TRIEE
!
!
    integer :: i, i1, i2, iageom, idlin1, idlin2, idlinv
    integer :: idlou1, idlou2, idlou3, idlou4, ier, iexcor, iret
    integer :: ino1, ino2, j, j1, j2, k
    integer :: nuno1, nuno2
!
    real(kind=8) :: d, dmin
    real(kind=8) :: mrot(3, 3), x1(3), x2(3)
!
    character(len=8) :: noma, m8blan
    character(len=8) :: nomno1, nomno2, nomo1, nomo2
    character(len=24) :: lisin1, lisin2, lisou1, lisou2
    character(len=24) :: valk(5)
    character(len=24) :: noeuma
    integer, pointer :: num_lisin1(:) => null()
    integer, pointer :: num_lisin2(:) => null()
!
! --- DEBUT
!
    call jemarq()
    lisin1 = lisi1z
    lisin2 = lisi2z
    lisou1 = liso1z
    lisou2 = liso2z
    noma = nomaz
    ier = 0
!
    noeuma = noma//'.NOMNOE'
    m8blan = '        '
    call jeveuo(noma//'.COORDO    .VALE', 'L', iageom)
!
! --- CONSTITUTION DE LA MATRICE DE ROTATION
!
    theta(1) = theta(1)*r8dgrd()
    theta(2) = theta(2)*r8dgrd()
    theta(3) = theta(3)*r8dgrd()
!
    call matrot(theta, mrot)
!
! --- CREATION SUR LA VOLATILE DES LISTES DE K8 LISOU1 ET LISOU2
! --- DE LONGUEUR LONLIS
!
    call jeexin(lisou1, iret)
    if (iret .ne. 0) then
        call jedetr(lisou1)
    end if
    call jeexin(lisou2, iret)
    if (iret .ne. 0) then
        call jedetr(lisou2)
    end if
    call wkvect(lisou1, 'V V K8', lonlis, idlou1)
    call wkvect(lisou2, 'V V K8', lonlis, idlou2)
!
    call jeveuo(lisin1, 'L', idlin1)
    call jeveuo(lisin2, 'L', idlin2)
!
! --- VECTEURS DE TRAVAIL
!
    call wkvect('&&PACOAP.LISOU3', 'V V K8', lonlis, idlou3)
    call wkvect('&&PACOAP.LISOU4', 'V V K8', lonlis, idlou4)
    call wkvect('&&PACOAP.LISINV', 'V V K8', lonlis, idlinv)
!
!     -- ON FABRIQUE UN OBJET QUI CONTIENDRA LES NUMEROS
!     -- DES NOEUDS DE LISIN1 ET LISIN2 :
    AS_ALLOCATE(vi=num_lisin1, size=lonlis)
    AS_ALLOCATE(vi=num_lisin2, size=lonlis)
    do k = 1, lonlis
        call jenonu(jexnom(noeuma, zk8(idlin1-1+k)), num_lisin1(k))
        call jenonu(jexnom(noeuma, zk8(idlin2-1+k)), num_lisin2(k))
    end do
!
! --- CONSTITUTION DE LA PREMIERE CORRESPONDANCE ENTRE LES LISTES
! --- DE NOEUDS LISIN1 ET LISIN2 ENTRE NO1 DONNE ET NO2 SELON LE
! --- CRITERE : NO2 = NO DANS LISIN2 / D(NO1,NO2) = MIN D(NO1,NO)
!
    do i1 = 1, lonlis
        nomno1 = zk8(idlin1+i1-1)
!       CALL JENONU(JEXNOM(NOEUMA,NOMNO1),NUNO1)
        nuno1 = num_lisin1(i1)
        call parotr(noma, iageom, nuno1, 0, centre, &
                    mrot, t, x1)
        dmin = r8gaem()
        j2 = 0
        do i2 = 1, lonlis
            nomo2 = zk8(idlin2+i2-1)
!         CALL JENONU(JEXNOM(NOEUMA,NOMO2),INO2)
            ino2 = num_lisin2(i2)
!         CALL PACOOR(NOMA,INO2,0,X2)
            x2(1) = zr(iageom-1+3*(ino2-1)+1)
            x2(2) = zr(iageom-1+3*(ino2-1)+2)
            x2(3) = zr(iageom-1+3*(ino2-1)+3)
            d = padist(3, x1, x2)
            if (d .lt. dmin) then
                dmin = d
                nomno2 = nomo2
                nuno2 = ino2
                j2 = i2
            end if
        end do
!
        if (j2 .eq. 0) then
            call utmess('F', 'MODELISA6_3', sk=nomno1)
        end if
!
        if (zk8(idlinv+j2-1) .eq. m8blan) then
            zk8(idlou1+i1-1) = nomno1
            zk8(idlou2+i1-1) = nomno2
            zk8(idlinv+j2-1) = nomno1
        else
            ier = ier+1
            valk(1) = nomno2
            valk(2) = nomno1
            valk(3) = zk8(idlinv+j2-1)
            call utmess('E', 'MODELISA8_77', nk=3, valk=valk)
        end if
!
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA6_4')
    end if
!
    do i = 1, lonlis
        zk8(idlinv+i-1) = m8blan
    end do
!
! --- CONSTITUTION DE LA SECONDE CORRESPONDANCE ENTRE LES LISTES
! --- DE NOEUDS LISIN1 ET LISIN2 ENTRE NO2 DONNE ET NO1 SELON LE
! --- CRITERE : NO1 = NO DANS LISIN1 / D(NO1,NO2) = MIN D(NO,NO2)
! --- LA CORRESPONDANCE EST DEFINIE PAR LA CONSTITUTION DES LISTES
! --- LISOU3 ET LISOU4.
!
    do i2 = 1, lonlis
        nomno2 = zk8(idlin2+i2-1)
        nuno2 = num_lisin2(i2)
        x2(1) = zr(iageom-1+3*(nuno2-1)+1)
        x2(2) = zr(iageom-1+3*(nuno2-1)+2)
        x2(3) = zr(iageom-1+3*(nuno2-1)+3)
        dmin = r8gaem()
        j1 = 0
        do i1 = 1, lonlis
            nomo1 = zk8(idlin1+i1-1)
            ino1 = num_lisin1(i1)
            call parotr(noma, iageom, ino1, 0, centre, &
                        mrot, t, x1)
            d = padist(3, x1, x2)
            if (d .lt. dmin) then
                dmin = d
                nomno1 = nomo1
                nuno1 = ino1
                j1 = i1
            end if
        end do
!
        if (j1 .eq. 0) then
            call utmess('F', 'MODELISA6_3', sk=nomno2)
        end if
!
        if (zk8(idlinv+j1-1) .eq. m8blan) then
            zk8(idlou3+i2-1) = nomno1
            zk8(idlou4+i2-1) = nomno2
            zk8(idlinv+j1-1) = nomno2
        else
            ier = ier+1
            valk(1) = nomno1
            valk(2) = nomno2
            valk(3) = zk8(idlinv+j1-1)
            call utmess('E', 'MODELISA8_77', nk=3, valk=valk)
        end if
!
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA6_4')
    end if
!
! --- VERIFICATION DE LA COHERENCE DES COUPLES FORMES D'UNE PART
! --- PAR LISOU1 ET LISOU2 ET D'AUTRE-PART DES COUPLES 'INVERSES'
! --- FORMES PAR LISOU3 ET LISOU4
!
    do i = 1, lonlis
        iexcor = 0
        do j = 1, lonlis
            if (zk8(idlou1+i-1) .eq. zk8(idlou3+j-1)) then
                iexcor = 1
                if (zk8(idlou2+i-1) .ne. zk8(idlou4+j-1)) then
                    ier = ier+1
                    valk(1) = lisin1
                    valk(2) = lisin2
                    valk(3) = zk8(idlou1+i-1)
                    valk(4) = zk8(idlou2+i-1)
                    valk(5) = zk8(idlou4+j-1)
                    call utmess('E', 'MODELISA8_87', nk=5, valk=valk)
                end if
            end if
        end do
!
        if (iexcor .eq. 0) then
            ier = ier+1
            valk(1) = lisin1
            valk(2) = lisin2
            valk(3) = zk8(idlou1+i-1)
            valk(4) = ' '
            valk(5) = ' '
            call utmess('E', 'MODELISA8_88', nk=5, valk=valk)
        end if
!
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA6_4')
    end if
!
! --- MENAGE
!
    call jedetr('&&PACOAP.LISOU3')
    call jedetr('&&PACOAP.LISOU4')
    call jedetr('&&PACOAP.LISINV')
    AS_DEALLOCATE(vi=num_lisin1)
    AS_DEALLOCATE(vi=num_lisin2)
!
    call jedema()
end subroutine
