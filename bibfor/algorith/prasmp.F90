! --------------------------------------------------------------------
! Copyright (C) 1991 - 201662025 - EDF R&D - www.code-aster.org
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

subroutine prasmp(option, nugene, tminbl, nomprn, modgen, &
                  tmnobl, tmadbl, knombl, inumbl, ssmax)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 13/10/92
!-----------------------------------------------------------------------
!  BUT:      < PREPARATION ASSEMBLAGE MATRICE PROJETEE >
!
!  PREPARER L'ASSEMBLAGE POUR UN LIGREL CORRESPONDANT AUX MATRICES
!   PROJETEES DES SOUS-STRUCTURE
!   ON CONSIDERE POUR L'ASSEMBLAGE UN LISTE GENERALE DES BLOC
!   ELEMENTAIRES A ASSEMBLER DANS UNE MATRICE STOCKEE PROFIL BLOC
!   (EN GENERAL MATRICE PROJETEE=1BLOC,MATRICE DE LIAISON=NBLOCS)
!   NOEUD TARDIF = NOEUD FICTIF SUPPORTANT UNE SOUS-STRUCTURE
!
!   ON REMPLIT TMNOBL TMADBL KNOMBL INUMBL
!
!-----------------------------------------------------------------------
!
! NOM----- / /:
!
! OPTION   /I/: NOM K11 DE L'OPTION D'ASSEMBLAGE
! NUGENE   /I/: NOM K14 DE LA NUMEROTATION GENERALISEE
! NOMPRN   /I/: NOM K8 DU LIGREL COURANT A TRAITER
! TMINBL   /I/: NOM K24 DE LA FAMILLE NOMMEE AU NOM DES LIGRELS
!               ET DONNANT POUR CHAQUE NOEUD TARDIF DU LIGREL
!               LE NUMERO DE SON 1 BLOC DANS LA LISTE GENERALE ET
!               LE NOMBRE DE BLOC
! MODGEN   /I/: NOM K8 MODELE_GENERALISE AMONT
! TMNOBL   /I/: NOM K24 DE LA FAMILLE NUMEROTE DONNANT POUR CHAQUE
!               TERME D'UN BLOC ELEMENTAIRE LE NUMERO DU BLOC ASSEMBLE
!               D'ARRIVE
! TMADBL   /I/: NOM K24 DE LA FAMILLE NUMEROTE DONNANT POUR CHAQUE
!               TERME D'UN BLOC ELEMENTAIRE LE RANG D'ARRIVEE
!               DANS LE  BLOC ASSEMBLE
! KNOMBL   /M/: VECTEUR DES NOM K24 DES OBJETS OU FAMILLE CONTENANT
!               LES BLOCS ELEMENTAIRES
! INUMBL   /M/: VECTEUR NUMERO  BLOCS ELEMENTAIRE DANS LEUR FAMILLE OU 0
!               LES BLOCS ELEMENTAIRES
! SSMAX    /M/: MAXIMUM DE LA VALEURS ABSOLUE DES TERMES TRAITES
!
!
!
#include "jeveux.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/maxblc.h"
#include "asterfort/maxblo.h"
#include "asterfort/mgutdm.h"
#include "asterfort/utmess.h"
#include "asterfort/nueq_chck.h"
!
!
!
    character(len=8) :: modgen, nomprn, nommcl, kbid
    character(len=14) :: nugene
    character(len=19) :: prgene, stomor
    character(len=9) :: rigopt, masopt, amoopt
    character(len=11) :: option, ricopt
    character(len=24) :: tmadbl, tmnobl, tminbl, knombl(*)
    character(len=24) :: valk
    character(len=10) :: adnom
    integer(kind=8) :: inumbl(*)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iad, ibid, ibl1, ieqc, inuc, inul, iret, i_ligr_ss
    integer(kind=8) :: j, lc, ll, llors, llprs, i_ligr
    integer(kind=8) :: ltadbl, ltinbl, ltnobl, nbcol, nblig, nbsst
    integer(kind=8) :: ntail, ntprno, nusst, ntria
    real(kind=8) :: ssmax
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
!-----------------------------------------------------------------------
    data rigopt, ricopt, masopt, amoopt/'RIGI_GENE', 'RIGI_GENE_C',&
     &                                 'MASS_GENE', 'AMOR_GENE'/
!-----------------------------------------------------------------------
!
    call jemarq()
    if (nomprn .ne. '&SOUSSTR') then
        goto 999
    end if
!
!------------------RECUPERATION DU NOMBRE DE SOUS-STRUCTURE-------------
    prgene = nugene//'.NUME'
    call nueq_chck(prgene)
    stomor = nugene//'.SMOS'
    call jenonu(jexnom(prgene//'.LILI', '&SOUSSTR'), i_ligr_ss)
    call jelira(jexnum(prgene//'.PRNO', i_ligr_ss), 'LONMAX', nbsst)
    nbsst = nbsst/2
!
!--------------------RECUPERATION DES CARACTERISTIQUES BLOCS------------
!
!
!------------------CREATION DU NOM A CONCATENER-------------------------
!   POUR RECUPERER LE NOM DES MATRICES PROJETEES
!
    if ((option .eq. rigopt) .or. (option .eq. ricopt)) then
        adnom = '.MAEL_RAID'
    else if (option .eq. masopt) then
        adnom = '.MAEL_MASS'
    else if (option .eq. amoopt) then
        adnom = '.MAEL_AMOR'
    end if
!
!---------------------REMPLISSAGE DES OBJETS DE TRAVAIL-----------------
!
!
    call jeveuo(prgene//'.NUEQ', 'L', vi=nueq)
    call jeveuo(stomor//'.SMDI', 'L', vi=smdi)
!
    call jenonu(jexnom('&&ASSGEN.REP.NOM.PROF', nomprn), ibid)
    call jeveuo(jexnum(tminbl, ibid), 'L', ltinbl)
    call jelira(jexnum(tminbl, ibid), 'LONMAX', ntprno)
    ntprno = ntprno/2
!
    call jenonu(jexnom(prgene//'.LILI', nomprn), i_ligr)
    call jeveuo(jexnum(prgene//'.ORIG', i_ligr), 'L', llors)
    call jeveuo(jexnum(prgene//'.PRNO', i_ligr), 'L', llprs)
!
!   BOUCLE SUR LES ELEMENTS DU LIGREL
!
    do j = 1, ntprno
!
!  RECUPERATION NUMERO SOUS-STRUCTURE
        nusst = zi(llors+j-1)
! RECUPERATION DU PROFIL BLOC ELEMENTAIRE
        ibl1 = zi(ltinbl+(j-1)*2)
!
! POUR UNE MATRICE ELEMENTAIRE BLOC IL FAUDRA FAIRE UNE BOUCLE
! SUR LES BLOC A ASSEMBLER (DO IBLOC=1,NBBL)
!
! RECUPERATION DU NOM DU MACR_ELEM AMONT
        kbid = '   '
        call mgutdm(modgen, kbid, nusst, 'NOM_MACR_ELEM', ibid, &
                    nommcl)
        knombl(ibl1) = nommcl//adnom//'_VALE'
        call jelira(knombl(ibl1), 'NMAXOC', ntria)
!
! VERIFICATION DE L'EXISTENCE DE LA MATRICE D'AMORTISSEMENT ASSOCIEE
! AU MACRO-ELEMENT
!
        if (option .eq. amoopt) then
            call jeexin(knombl(ibl1), iret)
            if (iret .eq. 0) then
                valk = nommcl
                call utmess('F', 'ALGORITH13_99', sk=valk)
            end if
        end if
!
        inumbl(ibl1) = 0
! TYPE DE LA MATRICE DU MACRO ELEMENT (IE REELLE OU COMPLEXE)
        if (option .eq. ricopt) then
            call maxblc(jexnum(nommcl//adnom//'_VALE           ', 1), ssmax)
        else
            call maxblo(jexnum(nommcl//adnom//'_VALE           ', 1), ssmax)
        end if
!  RECUPERATION DIMENSIONS ET NUMERO PREMIERE EQUATION DANS NUEQ
        nblig = zi(llprs+(j-1)*2+1)
        nbcol = nblig
        inul = zi(llprs+(j-1)*2)
        inuc = inul
!
! TAILLE BLOC MATRICE PROJETEE (STOCKAGE TRIANGLE SUPERIEUR)
!
        ntail = nbcol*(nbcol+1)/2
!
        call jecroc(jexnum(tmnobl, ibl1))
        call jeecra(jexnum(tmnobl, ibl1), 'LONMAX', ntail)
        call jeveuo(jexnum(tmnobl, ibl1), 'E', ltnobl)
        call jecroc(jexnum(tmadbl, ibl1))
        call jeecra(jexnum(tmadbl, ibl1), 'LONMAX', ntail)
        call jeveuo(jexnum(tmadbl, ibl1), 'E', ltadbl)
!
!     BOUCLES SUR LES TERMES DU BLOC ELEMENTAIRE
!    (TRIANGLE SUPERIEUR SEULEMENT)
        do ll = 1, nblig
            do lc = ll, nbcol
!    ADRESSE DANS BLOC ELEMENTAIRE
                iad = (lc-1)*lc/2+ll
!    NUMERO D'EQUATION DU TERME COURANT
                ieqc = nueq(1+(inuc-1)+(lc-1))
                zi(ltnobl+iad-1) = 1
                zi(ltadbl+iad-1) = smdi(ieqc)-(lc-ll)
!
            end do
        end do
        call jelibe(jexnum(tmnobl, ibl1))
        call jelibe(jexnum(tmadbl, ibl1))
    end do
!
999 continue
    call jedema()
end subroutine
