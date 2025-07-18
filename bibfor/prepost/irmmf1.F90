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
! person_in_charge: nicolas.sellenet at edf.fr
!
subroutine irmmf1(fid, nomamd, typent, nbrent, nbgrou, &
                  nomgen, nufaen, nomast, prefix, typgeo, &
                  nomtyp, nmatyp, infmed, ifm, nosdfu)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/irmmf2.h"
#include "asterfort/irmmf3.h"
#include "asterfort/jedetr.h"
#include "asterfort/wkvect.h"
!
    med_idt :: fid
    integer(kind=8) :: typent, nbrent, nbgrou
    integer(kind=8) :: nufaen(nbrent)
    integer(kind=8) :: typgeo(*), nmatyp(*)
    integer(kind=8) :: infmed
    integer(kind=8) :: ifm
    character(len=6) :: prefix
    character(len=8) :: nomast
    character(len=24) :: nomgen(*)
    character(len=8) :: nomtyp(*)
    character(len=*) :: nomamd
    character(len=8) :: nosdfu
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE - FORMAT MED - LES FAMILLES - 1
!
! --------------------------------------------------------------------------------------------------
!
!     L'ENSEMBLE DES FAMILLES EST L'INTERSECTION DE L'ENSEMBLE
!     DES GROUPES : UN ENTITE/MAILLE APPARAIT AU PLUS DANS 1 FAMILLE
!     TABLE  NUMEROS DES FAMILLES POUR LES ENTITES  <-> TABLE  DES COO
!     TABLES NUMEROS DES FAMILLES POUR MAILLE/TYPE <-> TABLES DES CNX
!     PAR CONVENTION, LES FAMILLES DE ENTITES SONT NUMEROTEES >0 ET LES
!     FAMILLES DE MAILLES SONT NUMEROTEES <0. LA FAMILLE NULLE EST
!     DESTINEE AUX ENTITES / ELEMENTS SANS FAMILLE.
!     ENTREE:
!   FID    : IDENTIFIANT DU FICHIER MED
!   NOMAMD : NOM DU MAILLAGE MED
!   TYPENT : TYPE D'ENTITES : 0, POUR DES NOEUDS, 1 POUR DES MAILLES
!   NBRENT : NOMBRE D'ENTITES A TRAITER
!   NBGROU : NOMBRE DE GROUPES D'ENTITES
!   NOMGEN : VECTEUR NOMS DES GROUPES D'ENTITES
!   NOMAST : NOM DU MAILLAGE ASTER
!   PREFIX : PREFIXE POUR LES TABLEAUX DES RENUMEROTATIONS
!   TYPGEO : TYPE MED POUR CHAQUE MAILLE
!   NOMTYP : NOM DES TYPES POUR CHAQUE MAILLE
!   NMATYP : NOMBRE DE MAILLES PAR TYPE
! TABLEAUX DE TRAVAIL
!   NUFAEN : NUMERO DE FAMILLE POUR CHAQUE ENTITE
!            PAR DEFAUT, L'ALLOCATION AVEC JEVEUX A TOUT MIS A 0. CELA
!            SIGNIFIE QUE LES ENTITES APPARTIENNENT A LA FAMILLE NULLE.
! DIVERS
!   INFMED : NIVEAU DES INFORMATIONS SPECIFIQUES A MED A IMPRIMER
!   IFM    : UNITE LOGIQUE DU FICHIER DE MESSAGE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRMMF1'
    integer(kind=8) :: tygeno
    integer(kind=8) :: iaux
    integer(kind=8) :: nbec
    integer(kind=8) :: adtabx, adnufa, adnogr, adnofe
    real(kind=8) :: raux
    character(len=7) :: saux07
    character(len=24) :: nufacr, nogrfa, nofaex, tabaux
!
! --------------------------------------------------------------------------------------------------
!
!
!====
! 1. S'IL EXISTE DES GROUPES, ON ALLOUE DES TABLEAUX DE TRAVAIL POUR
!    POUVOIR CONSTRUIRE LES FAMILLES
!====
!
    tygeno = 0
    if (infmed .ge. 2) then
!
        if (typent .eq. tygeno) then
            saux07 = 'NOEUDS '
        else
            saux07 = 'MAILLES'
        end if
        write (ifm, 100) saux07, saux07, nbrent, nbgrou
100     format(/, 'CONSTRUCTION DES FAMILLES DE ', a, /, '. NOMBRE DE ', a, ' :', i12, &
                /, '. NOMBRE DE GROUPES :', i5)
    end if
!
    if (nbgrou .ne. 0) then
!
!====
! 2. S'IL EXISTE DES GROUPES, ON ALLOUE DES TABLEAUX DE TRAVAIL POUR
!    POUVOIR CONSTRUIRE LES FAMILLES
!====
!                 12   345678   9012345678901234
        nufacr = '&&'//nompro//'.NU_FAM_CRE     '
        nogrfa = '&&'//nompro//'.NOM_GR_FAM     '
        nofaex = '&&'//nompro//'.NOM_FAM_EX     '
        tabaux = '&&'//nompro//'.TABL_AUXIL     '
!
!       VECTEUR NUMEROS DES FAMILLES D'ENTITES CREES = NBRENT MAX
        call wkvect(nufacr, 'V V I', nbrent, adnufa)
!
!       VECTEUR NOMS DES GROUPES D'ENTITES / FAMILLE = NB GRP MAX
        call wkvect(nogrfa, 'V V K80', nbgrou, adnogr)
!
!       VECTEUR NOMS DES FAMILLES = NB FAMILLE MAX
!       AU PIRE, IL Y A UNE FAMILLE PAR ENTITE. MAIS EN FAIT, C'EST
!       UNE PARTITION SELON LES GROUPES : IL Y EN A AU PIRE 2**NBGROU-1
!       ON CHOISIT DONC LE MIN DES 2
!       POURQUOI 2**NBGROU-1 ?
!       SOIT L'ENTITE APPARTIENT A 1 GROUPE  ==> NBGROU CHOIX
!       SOIT L'ENTITE APPARTIENT A 2 GROUPES ==> ARR(NBGROU,2) CHOIX
!       SOIT L'ENTITE APPARTIENT A 3 GROUPES ==> ARR(NBGROU,3) CHOIX
!       ...
!       SOIT L'ENTITE APPARTIENT A (NBGROU-1) GROUPES
!                                        ==> ARR(NBGROU,NBGROU-1) CHOIX
!       SOIT L'ENTITE APPARTIENT AUX NBGROU GROUPES ==> 1 CHOIX
!       AU TOTAL : NBGROU + ARR(NBGROU,2) + ARR(NBGROU,3) + ...
!                     ... + ARR(NBGROU,NBGROU-1) + 1
!       ON REMARQUE QUE CE SONT LES COEFFICIENTS D'UNE LIGNE DU TRIANGLE
!       DE PASCAL (HONNEUR AUX AUVERGNATS)
!       DONC LA SOMME VAUT (1+1)*NBGROU-1 - 1
!
        raux = log(dble(nbrent))/log(2.d0)
        if (nbgrou .lt. int(raux)) then
            iaux = 2**nbgrou-1
        else
            iaux = nbrent
        end if
        call wkvect(nofaex, 'V V K80', iaux, adnofe)
!
!       ON UTILISE DES TABLEAUX DE BITS POUR ENREGISTRER LA PRESENCE
!       D'UNE ENTITE DANS UN GROUPE : 30 PAR ENTIER INT*4 NBEC ENTIERS
!
        nbec = (nbgrou-1)/30+1
!       COLOSSAL ALLOC DU TABLEAU CROISE DES ENTITES X GROUPES
!
        call wkvect(tabaux, 'V V I', nbrent*nbec, adtabx)
!
!====
! 3. CREATION ET ECRITURE DES FAMILLES
!====
!
        if (nosdfu .eq. ' ') then
            call irmmf2(fid, nomamd, typent, nbrent, nbgrou, &
                        nomgen, nbec, nomast, prefix, typgeo, &
                        nomtyp, nmatyp, nufaen, zi(adnufa), zk80(adnogr), &
                        zk80(adnofe), zi(adtabx), infmed, ifm)
        else
            call irmmf3(fid, nomamd, typent, nbrent, nbgrou, &
                        nomgen, nbec, nomast, prefix, typgeo, &
                        nomtyp, nmatyp, nufaen, zi(adnufa), zk80(adnogr), &
                        zk80(adnofe), zi(adtabx), infmed, ifm, nosdfu)
        end if
!
        call jedetr(nufacr)
        call jedetr(nogrfa)
        call jedetr(nofaex)
        call jedetr(tabaux)
!
    end if
!
end subroutine
