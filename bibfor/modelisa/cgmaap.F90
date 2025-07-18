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

subroutine cgmaap(mofaz, iocc, nomaz, lismaz, nbma)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cgmaal.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: iocc, nbma
    character(len=*) :: mofaz, nomaz, lismaz
!
!       CGMAAP -- TRAITEMENT DE L'OPTION 'APPUI'
!                 DU MOT FACTEUR CREA_GROUP_MA DE
!                 LA COMMANDE DEFI_GROUP
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_MA'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISMAZ        - JXVAR - K24  - : NOM DE LA LISTE DE MAILLES OBTENUES
!                                   SUIVANT LE TYPE D'APPUI
!                                   (VOIR CI-DESSOUS)
!  NBMA          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!    TYPE D'APPUI:
!      - 'AU_MOINS_UN': UNE MAILLE EST RETENUE SI L'UN DE SES NOEUDS
!              EST DANS LA LISTE DES NOEUDS FOURNIS.
!      - 'TOUT': UNE MAILLE EST RETENUE SI TOUS SES NOEUDS
!              SONT DANS LA LISTE DES NOEUDS FOURNIS.
!      - 'SOMMET': UNE MAILLE EST RETENUE SI TOUS SES NOEUDS SOMMETS
!              SONT DANS LA LISTE DES NOEUDS FOURNIS.
!      - 'MAJORITE': UNE MAILLE EST RETENUE SI PLUS DE LA MOITIE DE
!             SES NOEUDS SONT DANS LA LISTE DES NOEUDS FOURNIS.
!
!
    integer(kind=8) :: nbmala, i, j, jmala, jco, iacnx, ilcnx, ii, nbmc
    integer(kind=8) ::  nbno, nno, jlmas, idlist, ima, n1
    integer(kind=8) ::  jnnma, nuno, j1, k, nbnot
    character(len=8) :: noma, motcle(2), tymocl(2), typma
    character(len=16) :: motfac, typapp
    character(len=24) :: listrv, lismai, mesnoe, lismas, lnnma
    integer(kind=8), pointer :: noeuds(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ================
    nbmala = 0
    motfac = mofaz
    noma = nomaz
    lismai = lismaz
    listrv = '&&CGMAAS.MAILLES_LACHE'
    mesnoe = '&&CGMAAL.NOEUDS'
    lismas = '&&CGMAAS.MAILLES_STRICT'
    lnnma = '&&CGMAAL.NBNO_MA'
!
! --  ON RECUPERE LES MAILLES DU TYPE D'APPUI 'AU_MOINS_UN'
!     (COMMMUN A TOUT TYPE D'APPUI)
    call cgmaal(motfac, iocc, noma, listrv, nbmala)
!
! --  RECUPERATION DU TYPE D'APPUI :
    call getvtx('CREA_GROUP_MA', 'TYPE_APPUI', iocc=iocc, scal=typapp, nbret=n1)
    call jeveuo(listrv, 'L', jmala)
!
! --- TYPE D'APPUI = 'AU_MOINS_UN'
!     ============================
    if (typapp .eq. 'AU_MOINS_UN') then
        nbma = nbmala
        jlmas = jmala
        goto 333
    else
        call wkvect(lismas, 'V V I', nbmala, jlmas)
    end if
!
! --- TYPE D'APPUI = 'TOUT' ET 'SOMMET'
!     ======================================
!
! --  RECUPERATIONS DES NOEUDS FOURNIS PAR L'UTILISATEUR
!     --------------------------------------------------
    nbmc = 2
    motcle(1) = 'NOEUD'
    tymocl(1) = 'NOEUD'
    motcle(2) = 'GROUP_NO'
    tymocl(2) = 'GROUP_NO'
    call reliem(' ', noma, 'NU_NOEUD', motfac, iocc, &
                nbmc, motcle, tymocl, mesnoe, nbno)
    call jeveuo(mesnoe, 'L', j1)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnot)
    AS_ALLOCATE(vi=noeuds, size=nbnot)
    do k = 1, nbno
        nuno = zi(j1-1+k)
        noeuds(nuno) = 1
    end do
!
    call jeveuo(noma//'.CONNEX', 'L', iacnx)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', ilcnx)
    call wkvect(lnnma, 'V V I', nbmala, jnnma)
!
!
    if (typapp .eq. 'SOMMET') then
!     --  ON DETERMINE LE NOMBRE DE NOEUDS 'SOMMET' POUR CHAQUE MAILLE
        call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
        do i = 1, nbmala
            ima = zi(jmala+i-1)
            call jenuno(jexnum('&CATA.TM.NOMTM', typmail(ima)), typma)
!
            if (typma(1:3) .eq. 'POI') then
                zi(jnnma+i-1) = 1
            else if (typma(1:3) .eq. 'SEG') then
                zi(jnnma+i-1) = 2
            else if (typma(1:4) .eq. 'QUAD') then
                zi(jnnma+i-1) = 4
            else if (typma(1:4) .eq. 'TRIA') then
                zi(jnnma+i-1) = 3
            else if (typma(1:5) .eq. 'TETRA') then
                zi(jnnma+i-1) = 4
            else if (typma(1:5) .eq. 'PENTA') then
                zi(jnnma+i-1) = 6
            else if (typma(1:5) .eq. 'PYRAM') then
                zi(jnnma+i-1) = 5
            else if (typma(1:4) .eq. 'HEXA') then
                zi(jnnma+i-1) = 8
            end if
        end do
!
    else if (typapp .eq. 'TOUT' .or. typapp .eq. 'MAJORITE') then
        do i = 1, nbmala
            ima = zi(jmala+i-1)
            zi(jnnma+i-1) = zi(ilcnx+ima)-zi(ilcnx+ima-1)
        end do
    else
        ASSERT(.false.)
    end if
!
! --  ON FILTRE LES MAILLES SUIVANT LES CRITERES RELATIFS
!     A CHAQUE TYPE D'APPUI:
!     ---------------------------------------------------
    nbma = 0
    do i = 1, nbmala
        nno = zi(jnnma+i-1)
        if (nno .eq. 0) goto 30
!
        ima = zi(jmala+i-1)
        jco = iacnx+zi(ilcnx+ima-1)-1
        ii = 0
        do j = 1, nno
            nuno = zi(jco+j-1)
            if (noeuds(nuno) .eq. 1) ii = ii+1
        end do
!
        if (typapp .eq. 'TOUT' .or. typapp .eq. 'SOMMET') then
            if (ii .eq. nno) then
                nbma = nbma+1
                zi(jlmas+nbma-1) = ima
            end if
        else if (typapp .eq. 'MAJORITE') then
            if (ii .ge. (nno+1)/2) then
                nbma = nbma+1
                zi(jlmas+nbma-1) = ima
            end if
        end if
30      continue
    end do
!
! --- CREATION ET REMPLISSAGE DU VECTEUR DE SORTIE
!     --------------------------------------------
333 continue
    if (nbma .gt. 0) then
        call wkvect(lismai, 'V V I', nbma, idlist)
        do i = 1, nbma
            zi(idlist+i-1) = zi(jlmas+i-1)
        end do
    end if
!
! --- FIN
!     ===
    call jedetr(listrv)
    call jedetr(lismas)
    call jedetr(mesnoe)
    call jedetr(lnnma)
    AS_DEALLOCATE(vi=noeuds)
!
    call jedema()
!
end subroutine
