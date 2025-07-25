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
subroutine cncinv(mail, lima, nlima, base, nomz)
!
! person_in_charge: jacques.pellet at edf.fr
!
!**********************************************************************
!   OPERATION REALISEE
!   ------------------
!     CONSTRUCTION DE LA TABLE DE CONNECTIVITE INVERSE, I.E. :
!
!     NOEUD --> LISTE DES MAILLES CONTENANT LE NOEUD
!
!   ATTENTION :
!   -----------
!      - SI LIMA EST FOURNI (NBMA > 0) , LA CONNECTIVITE INVERSE
!        RETOURNEE FAIT REFERENCE AUX INDICES DES MAILLES DANS LIMA:
!        NUMA=LIMA(CONINV(NUNO))
!      - SI LIMA N'EST PAS FOURNI (NBMA = 0) , LA CONNECTIVITE INVERSE
!        RETOURNEE FAIT REFERENCE AUX NUMEROS DES MAILLES DU MAILLAGE:
!        NUMA=CONINV(NUNO)
!      - CONVENTION : SI UN NOEUD NUNO EST ORPHELIN :
!        LONG(CONINV(NUNO))=1 ET CONINV(NUNO)(1)=0
!
!   ARGUMENTS EN ENTREE
!   ------------------
!     MAIL   : NOM DU MAILLAGE
!     LIMA   : LISTE DES NUMEROS DE MAILLES
!     NLIMA  : NOMBRE DE MAILLES DANS LIMA
!              SI NLIMA=0 ON PREND TOUT LE MAILLAGE
!     BASE   : BASE DE CREATION POUR NOMZ
!     NOMZ   : NOM DE L' OJB A CREER
!
!   ORGANISATION DE L'OBJET CREE DE NOM NOMZ :
!   --------------------------------------------
!     TYPE : XC V I ACCES(NUMEROTE) LONG(VARIABLE)
!
!**********************************************************************
!
    implicit none
!
!
!  FONCTIONS EXTERNES
!  ------------------
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
!
!  -----------------------------------------
!
!
!  ---------------------------------
!
! --- VARIABLES
    character(len=*) :: nomz
    character(len=24) :: nom
    character(len=8) :: mail
    character(len=1) :: base
    integer(kind=8) :: lima(*), nlima, nno, nma, ima, nare
    integer(kind=8) :: i, j, n, p0, p1, p2, p3, q1, q2, q3
    integer(kind=8), pointer :: indice(:) => null()
!
! --- LECTURE DONNEES
!
    call jemarq()
!
    nom = nomz
!
    call jeveuo(mail//'.DIME', 'L', p0)
    nno = zi(p0)
    if (nlima .eq. 0) then
        nma = zi(p0+2)
    else
        nma = nlima
    end if
!
    call jeveuo(mail//'.CONNEX', 'L', p1)
    call jeveuo(jexatr(mail//'.CONNEX', 'LONCUM'), 'L', p2)
!
    if ((nno .le. 0) .or. (nma .le. 0)) goto 90
!
! --- ALLOCATION OBJETS TEMPORAIRES
!
    AS_ALLOCATE(vi=indice, size=nma)
    call wkvect('&&CNCINV.NMAILLE', 'V V I', nno, q1)
    call wkvect('&&CNCINV.POINTEUR', 'V V I', nno+1, q2)
!
    if (nlima .eq. 0) then
!
        do i = 1, nma
            indice(i) = i
        end do
!
    else
!
        do i = 1, nma
            indice(i) = lima(i)
        end do
!
    end if
!
    do i = 1, nno
        zi(q1-1+i) = 0
    end do
!
! --- NOMBRE DE MAILLES POUR CHAQUE NOEUD
!
    do i = 1, nma
!
        ima = indice(i)
        p0 = zi(p2-1+ima)
        n = zi(p2+ima)-p0
        p0 = p1+p0-1
!
        do j = 1, n
            p3 = q1-1+zi(p0)
            p0 = p0+1
            zi(p3) = zi(p3)+1
        end do
    end do
!
! --- NOMBRE TOTAL D'ARETES NOEUD/MAILLE
!
    nare = 0
    zi(q2) = 0
!
    do i = 1, nno
!
        n = zi(q1-1+i)
        if (n .eq. 0) n = 1
!
        nare = nare+n
        zi(q2+i) = nare
!
    end do
!
! --- ALLOCATION DU GRAPHE NOEUD/MAILLE
!
    call jecrec(nom, base//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                nno)
    call jeecra(nom, 'LONT', nare)
!
    do i = 1, nno
!
        n = zi(q1-1+i)
        if (n .eq. 0) n = 1
        call jecroc(jexnum(nom, i))
        call jeecra(jexnum(nom, i), 'LONMAX', n)
!
    end do
!
! --- GRAPHE NOEUD/MAILLE
!
    call jeveuo(nom, 'E', q3)
!
    do i = 1, nare
        zi(q3-1+i) = 0
    end do
!
    do i = 1, nma
!
        ima = indice(i)
        p0 = zi(p2-1+ima)
        n = zi(p2+ima)-p0
        p0 = p1+p0-1
!
        do j = 1, n
!
            q1 = q2-1+zi(p0)
            p0 = p0+1
            zi(q3+zi(q1)) = i
            zi(q1) = zi(q1)+1
!
        end do
    end do
!
! --- DESALLOCATION
!
    AS_DEALLOCATE(vi=indice)
    call jedetr('&&CNCINV.NMAILLE')
    call jedetr('&&CNCINV.POINTEUR')
!
90  continue
!
    call jedema()
!
end subroutine
