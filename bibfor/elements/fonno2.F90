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

subroutine fonno2(macofo, noma, nbmac, nbnoff, nbnose, &
                  nbmax, noeu, tablev)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: noma, noeu
    character(len=19) :: macofo
    integer(kind=8) :: nbmac, nbnoff, nbnose, nbmax, tablev(2)
!
!      PARMI LES MAILLES CONNECTEES AU SEGMENT DU FOND, FILTRAGE DES
!          MAILLES CONNECTEES A 1 LEVRE (CAD AYANT UNE FACE LIBRE)
!          -> REMPLISSAGE DE TABLEV
!       ----------------------------------------------------------------
!    ENTREES
!       MACOFO : VECTEUR DES MAILLES (PRINCIPALES) CONNECTEES AU SEGMENT
!                DU FOND DE FISSURE COURANT
!       NOMA   : NOM DU MAILLAGE
!       NBMAC  : NOMBRE DE MAILLES CONNECTEES AU SEGMENT DU FOND ET DE
!                DE DIMENSION NDIM
!       NBNOFF : NOMBRE DE NOEUD EN FOND DE FISSURE
!       NBNOSE : NOMBRE DE NOEUD PAR SEGMENT
!       NBMAX  : NOMBRE DE NOEUDS MAX COMMUNS A DEUX MAILLES CONNEXES
!       NOEU   : NOM DU NOEUD SOMMET COURANT
!    SORTIE
!       TABLEV : VECTEUR CONTNANT LES NUMEROS DES DEUX MAILLES
!                CONNECTEES AU NOEUD SOMMET COURANT ET AUX LEVRES
!
!
    integer(kind=8) :: jmaco, iatyma, jno1, jno2, typ11, typ22
    integer(kind=8) :: inp, inq, inr, ins, nbno1, nbno2
    integer(kind=8) :: comp2, comp3, comp4
    character(len=8) :: typ1, typ2
    character(len=9) :: valk(1)
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
!
!     RECUPERATION DE L'ADRESSE DES TYPFON DE MAILLES
    call jeveuo(noma//'.TYPMAIL', 'L', iatyma)
!
!     RECUPERATION DU VECTEUR DES MAILLES CONNECTEES AU SEGMENT DU FOND
    call jeveuo(macofo, 'L', jmaco)
!
    comp4 = 0
    do inp = 1, nbmac
        comp3 = 0
        call jeveuo(jexnum(noma//'.CONNEX', zi(jmaco-1+inp)), 'L', jno1)
!       NOMBRE DE NOEUDS LA MAILLE
        typ11 = iatyma-1+zi(jmaco-1+inp)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(typ11)), typ1)
        call dismoi('NBNO_TYPMAIL', typ1, 'TYPE_MAILLE', repi=nbno1)
!       POUR CHAQUE MAILLE VOISINE (NOEUD FOND COMMUN ET MEME
!       DIMENSION TOPO
        do inq = 1, nbmac
            comp2 = 0
            call jeveuo(jexnum(noma//'.CONNEX', zi(jmaco-1+inq)), 'L', jno2)
            typ22 = iatyma-1+zi(jmaco-1+inq)
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(typ22)), typ2)
            call dismoi('NBNO_TYPMAIL', typ2, 'TYPE_MAILLE', repi=nbno2)
!         ON COMPTE LE NOMBRE DE NOEUDS COMMUNS AFIN D'ISOLER
!         LES MAILLES DE BORD
            do inr = 1, nbno1
                do ins = 1, nbno2
                    if (zi(jno1-1+inr) .eq. zi(jno2-1+ins)) then
                        comp2 = comp2+1
                    end if
                end do
            end do
!         SI LES DEUX MAILLES ONT DES NOEUDS EN COMMUN EN DEHORS DU
!         FOND MAIS PAS TOUS
            if (((nbnoff .eq. 1) .and. (comp2 .ne. nbno1) .and. (comp2 .ge. nbnose)) .or. &
                ((nbnoff .gt. 1) .and. (comp2 .ne. nbno1) .and. (comp2 .ge. nbmax)) .or. &
                (nbmac .eq. 1)) then
                comp3 = comp3+1
            end if
        end do
!       ON GARDE LES MAILLES CONNECTEES QU'A 1 SEULE AUTRE MAILLE
        if (comp3 .eq. 1) then
            comp4 = comp4+1
            ASSERT(comp4 .le. 2)
            tablev(comp4) = zi(jmaco-1+inp)
        end if
    end do
!
!     SI AUCUNE MAILLE DE CE TYPE N'EST TROUVE
    if (comp4 .eq. 0) then
        valk(1) = noeu
        call utmess('F', 'RUPTURE0_31', sk=valk(1))
    end if
    call jedema()
end subroutine
