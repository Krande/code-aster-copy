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

subroutine fonext(noma, cnxinv, jbasno, inoext, inoseg, &
                  nbnoff, jborl, jdirol, jnvdir, iseg)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cengra.h"
#include "asterfort/confac.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/xextre.h"
#include "asterfort/xfabor.h"
#include "asterfort/xnorme.h"
!
    integer(kind=8) :: jbasno, inoext, inoseg, nbnoff, jborl, jdirol
    integer(kind=8) :: jnvdir, iseg
    character(len=8) :: noma
    character(len=19) :: cnxinv
!
! FONCTION REALISEE:
!
!     CALCUL DU VECTEUR DE DIRECTION DE PROPAGATION AUX EXTREMITES DU
!     FOND DE FISSURE
!
!     ENTREES:
!        NOMA   : NOM DU MAILLAGE
!        CNXINV : CONNECTIVITE INVERSE
!        INOEXT : INDICE DU NOEUD EXTREMITE DU FOND
!        INOSEG : INDICE DU NOEUD APPARTENANT AU MEME SEGMENT DU FOND
!                 QUE INOEXT
!        NBNOFF : NOMBRE DE NOEUDS AU FOND DE FISSURE
!        JBORL  : ADRESSE DU VECTEUR PERMETTANT DE SAVOIR SI LE VECTEUR
!                  DE DIRECTION DE PROPAGATION A DEJA ETE RECALCULE OU
!                  NON AUX POINTS EXTREMITES DU FOND (POUR SAVOIR SI ON
!                  DOIT REMPLACER LA VALEUR EXISTANTE OU LA LUI AJOUTER)
!        JDIROL : ADRESSE DES VECTEURS DIRECTIONS DE PROPAGATION
!                  INITIAUX (CAD SANS MODIFICATION DES VECTEURS AUX
!                  POINTS EXTREMITES DU FOND)
!        JNVDIR : ADRESSE DU VECTEUR CONTENANT 0 OU 1 AUX POINTS
!                  EXTREMITES DU FOND:
!                  0: LE PRODUIT SCALAIRE ENTRE LA NORMALE A LA FACE DE
!                     BORD ET LE VDIR INITIAL ESI INFERIEUR A 0
!                  1: LE PRODUIT SCALAIRE EST SUPERIEUR OU EGAL A 0
!        ISEG   : INDICE DU SEGMENT DU FOND
!
!     ENTREE/SORTIE:
!        JBASNO : BASE LOCALE AUX NOEUDS DU FOND
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, ifa, ima, ino, itypma
    integer(kind=8) :: jconx2, jcoor, jmanoe
    integer(kind=8) :: nbf, nbfacb, nbno, ndime, nmaext, nmanoe, nuno
    integer(kind=8) :: nunoa, nunob, nunoc, numpt
    integer(kind=8) :: ibid3(12, 3), inobor(2), fa(6, 8)
    real(kind=8) :: coorg(3), vectn(12), norme, vect(3), proj
    character(len=8) :: typma
    aster_logical :: fabord, nofac
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!     -----------------------------------------------------------------
!
    call jemarq()
    numpt = 1
    if (iseg .ne. 1) numpt = nbnoff
!
!     RECUPERATION DES INFORMATIONS RELATIVES AU MAILLAGE
!
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call jeveuo(noma//'.COORDO    .VALE', 'L', jcoor)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
!     RECUPERATION DES MAILLES CONTENANT LE NOEUD EXTREMITE DU FOND
    call jelira(jexnum(cnxinv, inoext), 'LONMAX', nmanoe)
    call jeveuo(jexnum(cnxinv, inoext), 'L', jmanoe)
!
!     BOUCLE SUR LE NOMBRE DE MAILLES CONNECTEES AU NOEUD EXTREMITE
    do ima = 1, nmanoe
        nmaext = zi(jmanoe-1+(ima-1)+1)
!
        itypma = typmail(nmaext)
        call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
!
        call dismoi('DIM_TOPO', typma, 'TYPE_MAILLE', repi=ndime)
!
!       ON NE PREND QUE LES MAILLES EN 3D
        if (ndime .ne. 3) goto 100
!       CALCUL DU CENTRE DE GRAVITE DE LA MAILLE
        call cengra(noma, nmaext, coorg)
        call confac(typma, ibid3, ibid, fa, nbf)
        nbfacb = 0
        inobor(1) = 0
        inobor(2) = 0
!       BOUCLE SUR LE NOMBRE DE FACES DE LA MAILLE
        do ifa = 1, nbf
            nbno = 4
            if (fa(ifa, 4) .eq. 0) nbno = 3
            nofac = .false.
!         BOUCLE SUR LE NOMBRE DE NOEUDS DE LA FACE
            do ino = 1, nbno
                nuno = connex(zi(jconx2+nmaext-1)+fa(ifa, ino)-1)
                if (nuno .eq. inoseg) goto 110
                if (nuno .eq. inoext) nofac = .true.
!         FIN BOUCLE SUR LES NOEUDS
            end do
!
            if (nofac) then
                nunoa = connex(zi(jconx2+nmaext-1)+fa(ifa, 1)-1)
                nunob = connex(zi(jconx2+nmaext-1)+fa(ifa, 2)-1)
                nunoc = connex(zi(jconx2+nmaext-1)+fa(ifa, 3)-1)
!
!           ON VERIFIE SI LA FACE COURANTE EST UNE FACE DE BORD
                call xfabor(noma, cnxinv, nunoa, nunob, nunoc, &
                            fabord)
                if (fabord) then
                    call xnorme(numpt, inobor, vectn, nbfacb, nunoa, &
                                nunob, nunoc, jcoor, coorg)
!             ON VERIFIE QUE LA NORMALE A LA FACE N'EST PAS
!             COLINEAIRE AU VECTEUR NORMAL AU PLAN DE FISSURE
                    vect(1) = vectn(1+3*(nbfacb-1))
                    vect(2) = vectn(2+3*(nbfacb-1))
                    vect(3) = vectn(3+3*(nbfacb-1))
                    call normev(vect, norme)
!
                    proj = vect(1)*zr(jbasno-1+6*(numpt-1)+1)+vect(2)* &
                           zr(jbasno-1+6*(numpt-1)+2)+vect(3)*zr(jbasno-1+6* &
                                                                 (numpt-1)+3)
!
                    if (abs(proj) .ge. 0.95d0) then
                        nbfacb = nbfacb-1
                        goto 110
                    end if
!
                end if
            end if
!       FIN BOUCLE SUR LES FACES
110         continue
        end do
        if (nbfacb .ne. 0) then
            call xextre(inobor, vectn, nbfacb, jbasno, jborl, &
                        jdirol, jnvdir)
        end if
!     FIN BOUCLE SUR LES MAILLES
100     continue
    end do
!
    call normev(zr(jbasno-1+6*(numpt-1)+4), norme)
!
    call jedema()
end subroutine
