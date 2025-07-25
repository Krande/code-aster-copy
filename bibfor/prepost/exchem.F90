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
subroutine exchem(modloc, tcmp, nbc, nbsp, tvale, &
                  valcmp, taberr)
    implicit none
!
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/iposdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/ncpact.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: modloc(*), tcmp(*), nbc, taberr(*), nbsp
    real(kind=8) :: tvale(*), valcmp(*)
!
!*********************************************************************
!
!     OPERATION REALISEE
!     ------------------
!
!       EXTRACTION DES VALEURS D' UN ENSEMBLE DE COMPOSANTES SUR LES
!       NOEUDS D' UNE MAILLE DANS UN CHAM_ELEM
!
!     ARGUMENTS EN ENTREES
!     --------------------
!
!       MODLOC : VECTEUR MODE LOCALE DU TYPE D' ELEMENT DEFINI
!                SUR LA MAILLE
!
!                  (1) --> CODE (ICI TJS 3 : CAS DU CHAM_ELEM)
!
!                  (2) --> NUMERO DE LA GRANDEUR
!
!                  (3) --> NBR DE SCALAIRES UTILISES POUR DECRIRE
!                          LE CHAM_ELEM DE LA GRANDEUR SUR LA MAILLE
!
!                  (4) --> ACCES AU NBR DE POINTS UTILISES POUR LA
!                          DESCRIPTION
!
!                  (5) --> DEBUT DES DESCRIPTEURS PAR ENTIERS CODES
!
!
!       TCMP   : TABLE DES NUMEROS DES CMP ACTIVES POUR L' EXTRACTION
!       NBC    : NBR DE CMP ACTIVES
!       NBSP   : NBR DE SOUS-POINTS
!       TAVLE  : TABLE DES VALEURS DU CHAM_ELEM SUR LA MAILLE
!
!     ARGUMENTS EN SORTIE
!     -------------------
!
!       VALCMP : TABLE DES VALEURS DES CMP EXTRAITES
!
!*********************************************************************
!
    integer(kind=8) :: gd, nbpt, nbnmai, acpact, nbec, adesgd
    integer(kind=8) :: nbrcpa, adrnd, asgtnd, aposcp, poscmp, i, j, k
!
!   FONCTIONS EXTERNES
!   ------------------
!
!
!   -------------------------
!
!
!======================================================================
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call jemarq()
    gd = modloc(2)
!
    call jeveuo(jexnum('&CATA.GD.DESCRIGD', gd), 'L', adesgd)
!
    nbec = zi(adesgd+3-1)
    nbpt = modloc(4)
!
!     /* PAR HYP. D' APPEL : LE CHAMP EST REPRESENTE AUX NOEUDS */
!     /* DONC, MODLOC(4) < 0                                    */
!
    call ncpact(modloc(5), nbec, nbrcpa)
!
    if (nbpt .gt. 10000) then
!
!        /* CAS D' UNE REPRESENTATION VARIANT AVEC  LES NOEUDS */
!
        nbnmai = nbpt-10000
!
        call wkvect('&&EXCHEM.NBRCMPACTIVE', 'V V I', nbnmai, acpact)
        call wkvect('&&EXCHEM.ADRSGTNOEUD', 'V V I', nbnmai, asgtnd)
        call wkvect('&&EXCHEM.POSCMP', 'V V I', nbc*nbnmai, aposcp)
!
        zi(asgtnd+1-1) = 1
        zi(acpact+1-1) = nbrcpa
!
        do k = 1, nbc, 1
            i = 1
            zi(aposcp+k-1) = iposdg(modloc(5+(i-1)*nbec), tcmp(k))
!
        end do
!
        do i = 2, nbnmai, 1
!
            call ncpact(modloc(5+nbec*(i-1)), nbec, nbrcpa)
!
            zi(asgtnd+i-1) = zi(asgtnd+i-1-1)+nbrcpa*nbsp
            zi(acpact+i-1) = nbrcpa
!
            adrnd = (i-1)*nbc
!
            do k = 1, nbc, 1
!
                zi(aposcp+adrnd+k-1) = iposdg(modloc(5+(i-1)*nbec), tcmp( &
                                              k))
!
            end do
!
        end do
!
    else
!
!        /* CAS D' UNE REPRESENTATION CONSTANTES SUR LES NOEUDS */
!
        nbnmai = nbpt
!
        call wkvect('&&EXCHEM.NBRCMPACTIVE', 'V V I', nbnmai, acpact)
        call wkvect('&&EXCHEM.ADRSGTNOEUD', 'V V I', nbnmai, asgtnd)
        call wkvect('&&EXCHEM.POSCMP', 'V V I', nbc*nbnmai, aposcp)
!
        do i = 1, nbnmai, 1
!
            zi(asgtnd+i-1) = (i-1)*nbrcpa*nbsp+1
            zi(acpact+i-1) = nbrcpa
!
            adrnd = (i-1)*nbc
!
            do k = 1, nbc, 1
!
                zi(aposcp+adrnd+k-1) = iposdg(modloc(5), tcmp(k))
!
            end do
!
        end do
!
    end if
!
    do i = 1, nbnmai, 1
!
        adrnd = zi(asgtnd+i-1)
        nbrcpa = zi(acpact+i-1)
!
        do j = 1, nbsp, 1
!
            do k = 1, nbc, 1
!
                poscmp = zi(aposcp+(i-1)*nbc+k-1)
!
                if (poscmp .gt. 0) then
!
                    valcmp(((i-1)*nbsp+j-1)*nbc+k) = tvale(adrnd+(j-1)*nbrcpa+poscmp-1)
!
                    taberr(k) = 1
!
                else
!
                    valcmp(((i-1)*nbsp+j-1)*nbc+k) = r8vide()
!
                    taberr(k) = 0
!
                end if
!
            end do
!
        end do
!
    end do
!
    call jedetr('&&EXCHEM.NBRCMPACTIVE')
    call jedetr('&&EXCHEM.ADRSGTNOEUD')
    call jedetr('&&EXCHEM.POSCMP')
!
    call jedema()
end subroutine
