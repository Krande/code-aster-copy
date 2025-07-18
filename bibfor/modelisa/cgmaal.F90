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

subroutine cgmaal(mofaz, iocc, nomaz, lismaz, nbma)
    implicit none
#include "jeveux.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: iocc, nbma
    character(len=*) :: mofaz, nomaz, lismaz
!
!       CGMAAL -- TRAITEMENT DE L'OPTION APPUI_LACHE
!                 DU MOT FACTEUR CREA_GROUP_MA DE
!                 LA COMMANDE DEFI_GROUP
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_MA'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISMAZ        - JXVAR - K24  - : NOM DE LA LISTE DE MAILLES QUI
!                                   CONTIENNENT NOEUD UTILISATEUR
!  NBMA          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
    integer(kind=8) :: nbmc, nbno, nci, adrvlc, acncin, i, j, nbmat, ityp
    integer(kind=8) :: nuno, jadr, numa, idlist, jnoeu, idlima
    character(len=8) :: noma, motcle(2), tymocl(2)
    character(len=16) :: motfac
    character(len=24) :: mesnoe, lismai, listrv, ncncin
!     -----------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    motfac = mofaz
    noma = nomaz
    lismai = lismaz
    listrv = '&&CGMAAL.MAILLES_TRAV'
    mesnoe = '&&CGMAAL.NOEUDS'
!
! --- RECUPERATION DES NOEUDS :
!     -----------------------
    nbmc = 2
    motcle(1) = 'NOEUD'
    tymocl(1) = 'NOEUD'
    motcle(2) = 'GROUP_NO'
    tymocl(2) = 'GROUP_NO'
    call reliem(' ', noma, 'NU_NOEUD', motfac, iocc, &
                nbmc, motcle, tymocl, mesnoe, nbno)
    call jeveuo(mesnoe, 'L', jnoeu)
!
    ncncin = '&&OP0104.CONNECINVERSE  '
    call jeexin(ncncin, nci)
    if (nci .eq. 0) call cncinv(noma, [0], 0, 'V', ncncin)
!
    call jeveuo(jexatr(ncncin, 'LONCUM'), 'L', adrvlc)
    call jeveuo(jexnum(ncncin, 1), 'L', acncin)
!
!
! --- RECUPERATION DU NOMBRE DE MAILLES DU MAILLAGE :
!     ---------------------------------------------
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmat)
!
! --- RECUPERATION DU TYPE DES MAILLES DU MAILLAGE :
!     --------------------------------------------
    call jeveuo(noma//'.TYPMAIL', 'L', ityp)
!
! --- ALLOCATION D'UN VECTEUR DE TRAVAIL :
!     ----------------------------------
    call wkvect(listrv, 'V V I', nbmat, idlima)
!
! --- TRAITEMENT DES NOEUDS "UTILISATEUR" :
!     -----------------------------------
    do i = 1, nbno
        nuno = zi(jnoeu+i-1)
        nbma = zi(adrvlc+nuno+1-1)-zi(adrvlc+nuno-1)
        jadr = zi(adrvlc+nuno-1)
        do j = 1, nbma
            numa = zi(acncin+jadr-1+j-1)
            zi(idlima+numa-1) = 1
        end do
    end do
!
    nbma = 0
    do i = 1, nbmat
        if (zi(idlima+i-1) .eq. 1) nbma = nbma+1
    end do
!
! --- ALLOCATION DU VECTEUR DES NOMS DES MAILLES CONTENANT
!     LES NOEUDS LISTES :
!     -----------------
    if (nbma .gt. 0) then
        call wkvect(lismai, 'V V I', nbma, idlist)
        nbma = 0
        do i = 1, nbmat
            if (zi(idlima+i-1) .eq. 1) then
                nbma = nbma+1
                zi(idlist+nbma-1) = i
            end if
        end do
    end if
!
    call jedetr(listrv)
    call jedetr(mesnoe)
!
    call jedema()
!
end subroutine
