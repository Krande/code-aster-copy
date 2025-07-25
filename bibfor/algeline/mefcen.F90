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
subroutine mefcen(caelem, iequiv, nbcyl, nbz, irot, &
                  numnog, nbnog, nummag, numgrp, coor, &
                  cent, req, xint, yint, zint, &
                  rint, nbgrp)
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: iequiv, nbcyl, numnog(*), nbnog(*), nummag(*)
    integer(kind=8) :: numgrp(*), irot(3), nbgrp, nbz
    real(kind=8) :: xint(nbcyl), yint(nbcyl), zint(nbz, nbgrp), coor(*)
    real(kind=8) :: cent(2*nbcyl), req(nbgrp), rint(nbcyl)
    character(len=19) :: caelem
!     RECUPERATION DANS LES CARTES ELEMENTS DES RAYONS DANS LE CAS OU
!     IL N Y A PAS DE GROUPES D EQUIVALENCE ET DANS LES DONNEES DANS
!     LE CAS OU IL Y A DES GROUPES D EQUIVALENCE - RECUPERATION DES
!     COORDONNEES DES CENTRES DES CYLINDRES REELS DANS LES DONNEES
!     OPERATEUR APPELANT : OP0144 , FLUST3
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : CAELEM : NOM DU CONCEPT DE TYPE CARA_ELEM
! IN  : IEQUIV : INDICE D EXISTANCE DES GROUPES D EQUIVALENCE
! IN  : NBCYL  : NOMBRE DE CYLINDRES REELS
! IN  : NBZ    : NOMBRE DE NOEUDS DE LA DISCRETISATION AXIALE
! IN  : IROT   : INDICE DE PERMUTATION CIRCULAIRE DU CHANGEMENT DE
!                REPERE
! IN  : NUMNOG : TABLEAU DES ADRESSES DES NUMEROS DES NOEUDS DES
!                CYLINDRES
! IN  : NBNOG  : TABLEAU DU NOMBRE DE NOEUDS DE CHAQUE CYLINDRE
! IN  : NUMMAG : TABLEAU DES ADRESSES DES NUMEROS DES MAILLES DES
!                CYLINDRES
! IN  : NUMGRP : INDICES DES GROUPES D EQUIVALENCE
! IN  : COOR   : COORDONNEES DES NOEUDS DU MAILLAGE
! IN  : CENT   : COORDONNEES DES CENTRES DONNES DANS LA COMMANDE
! IN  : REQ    : RAYONS DES CYLINDRES DONNES DANS LA COMMANDE
! OUT : XINT   : COORDONNEES 'X' DANS LE REPERE AXIAL DES CENTRES DES
!                CYLINDRES
! OUT : YINT   : COORDONNEES 'Y' DANS LE REPERE AXIAL DES CENTRES DES
!                CYLINDRES
! OUT : ZINT   : COORDONNEES 'Z' DANS LE REPERE AXIAL DES NOEUDS DES
!                CYLINDRES
! OUT : RINT   : RAYONS DES CYLINDRES
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, iret, rangr1
    character(len=19) :: carte, carsd
    character(len=3) :: note
!     ------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iad, icesc, icesd, icesl, icmp
    integer(kind=8) :: npmax, numma, numno1, numno2
    real(kind=8) :: epsit
    real(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    epsit = 1.d-5
!
!
! --- COORDONNEES DES CENTRES DES CYLINDRES
! --- CAS OU IL N Y A PAS DE GROUPES D EQUIVALENCE
!
    if (iequiv .eq. 0) then
        do i = 1, nbcyl
            numno1 = zi(numnog(i))
            xint(i) = coor((numno1-1)*3+irot(1))
            yint(i) = coor((numno1-1)*3+irot(2))
            zint(1, numgrp(i)) = coor((numno1-1)*3+irot(3))
            do j = 2, nbnog(i)
                numno2 = zi(numnog(i)+j-1)
                if (abs(coor((numno1-1)*3+irot(1))-coor((numno2-1)*3+irot(1))) &
                    .gt. epsit .or. &
                    abs(coor((numno1-1)*3+irot(2))-coor((numno2-1)*3+irot(2))) &
                    .gt. epsit) then
                    write (note(1:3), '(I3.3)') i
                    call utmess('F', 'ALGELINE_73', sk=note)
                end if
                zint(j, numgrp(i)) = coor((numno2-1)*3+irot(3))
            end do
        end do
!
!
! --- COORDONNEES DES CENTRES DES CYLINDRES
! --- CAS OU IL Y A DES GROUPES D EQUIVALENCE
!
    else if (iequiv .eq. 1) then
        do i = 1, nbcyl
            xint(i) = cent(2*(i-1)+1)
            yint(i) = cent(2*(i-1)+2)
            do j = 1, nbnog(numgrp(i))
                numno2 = zi(numnog(numgrp(i))+j-1)
                zint(j, numgrp(i)) = coor((numno2-1)*3+irot(3))
            end do
        end do
    end if
!
!
! --- RAYONS DES CYLINDRES
! --- CAS OU IL Y A DES GROUPES D EQUIVALENCE
!
    if (iequiv .eq. 1) then
        do i = 1, nbcyl
            rint(i) = req(numgrp(i))
        end do
!
! --- RAYONS DES CYLINDRES
! --- CAS OU IL N Y A PAS DES GROUPES D EQUIVALENCE
!
    else if (iequiv .eq. 0) then
!CC ON RECUPERE LA CARTE ET ON LA TRANSFORME EN CHAMELEM_S
        carte = caelem(1:8)//'.CARGEOPO'
        carsd = '&&MEFCEN.CARGEOPO'
        call carces(carte, 'ELEM', ' ', 'G', carsd, &
                    'A', iret)
        ASSERT(iret .eq. 0)
!
! --- RECUPERATION DE LA GRANDEUR (ICI R1)  ---
! --- REFERENCEE PAR LA CARTE CARGEOPO           ---
!
        call jeveuo(carsd//'.CESC', 'L', icesc)
        call jeveuo(carsd//'.CESD', 'L', icesd)
        call jeveuo(carsd//'.CESL', 'L', icesl)
        call jeveuo(carsd//'.CESV', 'L', vr=cesv)
!
        call jeveuo(carte//'.DESC', 'L', vi=desc)
        call jeveuo(jexnom('&CATA.GD.NOMCMP', 'CAGEPO_R'), 'L', icmp)
!
        call jelira(jexnum('&CATA.GD.NOMCMP', desc(1)), 'LONMAX', npmax)
!
        rangr1 = indik8(zk8(icmp), 'R1      ', 1, npmax)
!
! ---    DEBUT DE LA BOUCLE SUR LES CYLINDRES
! ---    ON RECHERCHE LA MAILLE ASSOCIEE AU PREMIER NOEUDS DE CHAQUE
! ---    CYLINDRE, ET ON LIT LA CARTE ELEMENT QUI LUI CORRESPOND
        do i = 1, nbcyl
            numma = zi(nummag(i))
!
            call cesexi('C', icesd, icesl, numma, 1, &
                        1, rangr1, iad)
            if (iad .gt. 0) then
! ---       RECUPERATION DU RAYON DE LA PREMIERE MAILLE DE CHAQUE
! ---       CYLINDRE
!
                rint(i) = cesv(abs(iad))
!
            else
                call utmess('F', 'ALGELINE_75')
            end if
        end do
    end if
!
! --- MENAGE
    call detrsd('CHAM_ELEM_S', '&&MEFCEN.CARGEOPO')
!
    call jedema()
end subroutine
