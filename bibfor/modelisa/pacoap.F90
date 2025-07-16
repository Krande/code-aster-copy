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
#include "asterfort/jeveuo.h"
#include "asterfort/matrot.h"
#include "asterfort/parotr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: lonlis
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
    integer(kind=8) :: i1, i2, iageom, idlin1, idlin2
    integer(kind=8) :: idlou1, idlou2, ier, iret
    integer(kind=8) :: ino2, jmin, k
    integer(kind=8) :: nuno1, nuno2
!
    real(kind=8) :: dsquared, dmin
    real(kind=8) :: mrot(3, 3), x1(3), dx(3)
!
    character(len=8) :: noma, m8blan
    character(len=8) :: nomno1, nomno2, nomo2
    character(len=24) :: lisin1, lisin2, lisou1, lisou2
    character(len=24) :: valk(5)
    integer(kind=8), pointer :: num_lisin1(:) => null()
    integer(kind=8), pointer :: num_lisin2(:) => null()
    character(len=8), pointer :: lisinv(:) => null()

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
    AS_ALLOCATE(vk8=lisinv, size=lonlis)
    AS_ALLOCATE(vi=num_lisin1, size=lonlis)
    AS_ALLOCATE(vi=num_lisin2, size=lonlis)
!
!     -- ON FABRIQUE UN OBJET QUI CONTIENDRA LES NUMEROS
!     -- DES NOEUDS DE LISIN1 ET LISIN2 :
!
    do k = 1, lonlis
        num_lisin1(k) = char8_to_int(zk8(idlin1-1+k))
        num_lisin2(k) = char8_to_int(zk8(idlin2-1+k))
    end do
!
! --- CONSTITUTION DE LA CORRESPONDANCE ENTRE LES LISTES
! --- DE NOEUDS LISIN1 ET LISIN2 ENTRE NO1 DONNE ET NO2 SELON LE
! --- CRITERE : NO2 = NO DANS LISIN2 / D(NO1,NO2) = MIN D(NO1,NO)
    do i1 = 1, lonlis
        nomno1 = zk8(idlin1+i1-1)
        nuno1 = num_lisin1(i1)
        call parotr(noma, iageom, nuno1, 0, centre, mrot, t, x1)
        dmin = r8gaem()
        jmin = 0
        do i2 = 1, lonlis
            nomo2 = zk8(idlin2+i2-1)
            ino2 = num_lisin2(i2)
            dx(1) = x1(1)-zr(iageom-1+3*(ino2-1)+1)
            dx(2) = x1(2)-zr(iageom-1+3*(ino2-1)+2)
            dx(3) = x1(3)-zr(iageom-1+3*(ino2-1)+3)
            dsquared = dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
            if (dsquared .lt. dmin) then
                dmin = dsquared
                nomno2 = nomo2
                nuno2 = ino2
                jmin = i2
            end if
        end do
!
        if (jmin .eq. 0) then
            call utmess('F', 'MODELISA6_3', sk=nomno1)
        end if
!
        if (lisinv(jmin) .eq. m8blan) then
            zk8(idlou1+i1-1) = nomno1
            zk8(idlou2+i1-1) = nomno2
            lisinv(jmin) = nomno1
        else
            ier = ier+1
            valk(1) = nomno2
            valk(2) = nomno1
            valk(3) = lisinv(jmin)
            call utmess('E', 'MODELISA8_77', nk=3, valk=valk)
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
    AS_DEALLOCATE(vi=num_lisin1)
    AS_DEALLOCATE(vi=num_lisin2)
    AS_DEALLOCATE(vk8=lisinv)

!
    call jedema()
end subroutine
