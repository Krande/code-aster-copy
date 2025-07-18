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
subroutine xoriff(info, nfon, jfono, jbaso, jtailo, &
                  nmafon, listpt, goinop, jfon, jnofaf, &
                  jbas, jtail, fonmul, nbfond)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/padist.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xffcr.h"
#include "asterfort/xffext.h"
    integer(kind=8) :: nfon, jfono, jbaso, jtailo, nmafon, jfon, jbas, jtail
    integer(kind=8) :: jnofaf, nbfond
    character(len=19) :: info, listpt
    character(len=24) :: fonmul
    aster_logical :: goinop
!
! ----------------------------------------------------------------------
!       ORIENTATION DES POINTS DU FOND DE FISSURE DANS LE CADRE DE XFEM
!
!  ENTRESS :
!     INFO  :   NOM DU VECTEUR INFO DE LA SD
!     JFONO :   ADRESSE DES POINTS DU FOND DE FISSURE DÉSORDONNÉS
!     JBASO :   ADRESSE DES DIRECTIONS DE PROPAGATION DÉSORDONNÉES
!     JTAILO:   ADRESSE DES TAILLES MAXIMALES DE MAILLES DÉSORDONNÉES
!     NFON  :   NOMBRE DE POINTS DU FOND DE FISSURE DÉSORDONNÉS
!     LISTPT:   LISTE DES INDICES DES POINTS DU FOND DÉSORDONNÉS
!     GOINOP :  .TRUE.  SI  OPOO10 AVEC UPWIND-SIMPLEXE/GRILLE/3D
!               .FALSE. SINON
!
!  SORTIES :
!     JFON  :  ADRESSE DES POINTS DU FOND DE FISSURE ORDONNÉS
!     JNOFAF:  ADRESSE DES NUMERO DES NOEUDS DES FACES DES ELEMENTS
!              PARENTS QUI CONTIENNENT LES POINTS DU FOND DE FISSURE
!              ORDONNES
!     JBAS  :  ADRESSE DES DIRECTIONS DE PROPAGATION ORDONNÉES
!     JTAIL :  ADRESSE DES TAILLES MAXIMALES DE MAILLES ORDONNÉES
!     FONMUL:  VECTEUR CONTENANT LE DEBUT ET L'ARRIVEE DE
!              CHAQUE FOND DE FISSURE
!     NBFOND:  NOMBRE DE FONDS DE FISSURE
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: indice, indicm, indipt, ima, ipt, iptext
    integer(kind=8) :: jfonmu, jinfo, jlistp, jptext, jtabpt, k, nbptex
    real(kind=8) :: absc, m(3), p(3)
    character(len=19) :: ptextr, tabpt, typfon
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
!     RECHERCHE DES POINTS EXTREMITES DU FOND DE FISSURE
!
    call jeveuo(info, 'L', jinfo)
    ptextr = '&&XORIFF.PTEXTR'
    call xffext(jinfo, nfon, nmafon, listpt, ptextr, &
                nbptex)
    typfon = ' '
    if (zk16(jinfo-1+3) .eq. 'FERME') typfon = 'FERME'
!
    call jeveuo(listpt, 'L', jlistp)
    call jeveuo(ptextr, 'L', jptext)
    call jeveuo(fonmul, 'L', jfonmu)
!
!     VECTEUR DES INDICES DES POINTS DU FOND ORDONNES
    tabpt = '&&XORIFF.TABPT'
    call wkvect(tabpt, 'V V I', nfon, jtabpt)
!
!     INDICE DU PREMIER POINT
    indipt = zi(jptext-1+1)
    zi(jtabpt-1+1) = indipt
    zi(jptext-1+1) = 0
    zi(jfonmu-1+1) = 1
    zr(jfono-1+11*(1-1)+4) = 0.d0
!
    nbfond = 1
!
    do ipt = 1, nfon-1
        do ima = 1, nmafon
!
            if (zi(jlistp-1+2*(ima-1)+2) .eq. 0) goto 11
!
!---      LE PREMIER INDICE CORRESPOND A CELUI RECHERCHE
            if (zi(jlistp-1+2*(ima-1)+1) .eq. zi(jtabpt-1+ipt)) then
                if (ipt .gt. 2) then
                    if (zi(jtabpt-1+ipt-2) .eq. zi(jlistp-1+2*(ima-1)+2)) then
!               DOUBLON DANS LA LISTE
                        zi(jlistp-1+2*(ima-1)+2) = 0
                        goto 11
                    end if
                end if
!
                zi(jtabpt-1+ipt+1) = zi(jlistp-1+2*(ima-1)+2)
                zi(jlistp-1+2*(ima-1)+2) = 0
!
!           CALCUL DE L'ABSCISSE CURVILIGNE
                indipt = zi(jtabpt-1+ipt+1)
                indicm = zi(jtabpt-1+ipt)
                do k = 1, 3
                    p(k) = zr(jfono-1+11*(indipt-1)+k)
                    m(k) = zr(jfono-1+11*(indicm-1)+k)
                end do
                absc = zr(jfono-1+11*(indicm-1)+4)
                zr(jfono-1+11*(indipt-1)+4) = absc+padist(3, m, p)
!
!---      LE DEUXIEME INDICE CORRESPOND A CELUI RECHERCHE
            elseif (zi(jlistp-1+2*(ima-1)+2) .eq. zi(jtabpt-1+ipt)) &
                then
                if (ipt .gt. 2) then
                    if (zi(jtabpt-1+ipt-2) .eq. zi(jlistp-1+2*(ima-1)+1)) then
!               DOUBLON DANS LA LISTE
                        zi(jlistp-1+2*(ima-1)+2) = 0
                        goto 11
                    end if
                end if
!
                zi(jtabpt-1+ipt+1) = zi(jlistp-1+2*(ima-1)+1)
                zi(jlistp-1+2*(ima-1)+2) = 0
!
!           CALCUL DE L'ABSCISSE CURVILIGNE
                indipt = zi(jtabpt-1+ipt+1)
                indicm = zi(jtabpt-1+ipt)
                do k = 1, 3
                    p(k) = zr(jfono-1+11*(indipt-1)+k)
                    m(k) = zr(jfono-1+11*(indicm-1)+k)
                end do
                absc = zr(jfono-1+11*(indicm-1)+4)
                zr(jfono-1+11*(indipt-1)+4) = absc+padist(3, m, p)
!
            end if
11          continue
        end do
!
!       ON N'A PAS TROUVE DE POINT A ASSOCIER A IPT: C'EST UN POINT
!       EXTREMITE ( CAS DES FONDS MULTIPLES )
        if (zi(jtabpt-1+ipt+1) .eq. 0) then
!
!         PRESENCE DE PLUSIEURS FONDS FERMES INTERDIT
            if (typfon .eq. 'FERME') call utmess('F', 'XFEM_20')
!
            indice = 0
!
!         VERIFICATION QUE LE DERNIER POINT EST UNE EXTREMITE DU FOND
            do iptext = 1, nbptex
                if (zi(jptext-1+iptext) .eq. zi(jtabpt-1+ipt)) then
                    zi(jfonmu-1+2*(nbfond-1)+2) = ipt
                    indice = 1
                    goto 13
                end if
            end do
13          continue
!
            ASSERT(indice .ne. 0)
            zi(jptext-1+iptext) = 0
!
!         RECHERCHE D'UN NOUVEAU POINT D'EXTREMITE POUR DEBUTER LE
!         NOUVEAU FOND DE FISSURE
            do iptext = 1, nbptex
                if (zi(jptext-1+iptext) .ne. 0) then
!
                    indipt = zi(jptext-1+iptext)
                    zi(jtabpt-1+ipt+1) = indipt
!
                    zr(jfono-1+11*(indipt-1)+4) = 0.d0
                    zi(jptext-1+iptext) = 0
                    nbfond = nbfond+1
                    zi(jfonmu-1+2*(nbfond-1)+1) = ipt+1
!
                    goto 15
                end if
            end do
!
!         PRESENCE DE FONDS OUVERTS ET DE FONDS FERMES INTERDIT
            call utmess('F', 'XFEM_21')
!
        end if
15      continue
    end do
!
    zi(jfonmu-1+2*(nbfond-1)+2) = nfon
!
    if (typfon .eq. 'FERME') then
        zi(jfonmu-1+2*(nbfond-1)+2) = nfon+1
    end if
!
    call utmess('I', 'XFEM_34', si=nbfond)
!
!     ORDONNANCEMENT DE FONDFISS, DE NOFACPTFON, DE BASEFOND
!     ET DE FOND.TAILLE_R
    call xffcr(nfon, jfono, jbaso, jtailo, jtabpt, &
               typfon, jfon, jnofaf, jbas, jtail)
!
    if (goinop) then
        call jedetr(ptextr)
        call jedetr(tabpt)
    end if
    call jedema()
end subroutine
