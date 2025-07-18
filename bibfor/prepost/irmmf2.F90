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
subroutine irmmf2(fid, nomamd, typent, nbrent, nbgrou, &
                  nomgen, nbec, nomast, prefix, typgeo, &
                  nomtyp, nmatyp, nufaen, nufacr, nogrfa, &
                  nofaex, tabaux, infmed, ifm)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/as_mfacre.h"
#include "asterfort/as_mmhfnw.h"
#include "asterfort/desgfa.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mdnofa.h"
#include "asterfort/nomgfa.h"
#include "asterfort/setgfa.h"
#include "asterfort/utmess.h"
!
    med_idt :: fid
    integer(kind=8) :: typgeo(*), nmatyp(*)
    integer(kind=8) :: typent, nbrent, nbgrou
    integer(kind=8) :: nbec
    integer(kind=8) :: nufaen(nbrent), nufacr(nbrent), tabaux(*)
    integer(kind=8) :: infmed
    integer(kind=8) :: ifm
    character(len=6) :: prefix
    character(len=8) :: nomast
    character(len=24) :: nomgen(*)
    character(len=8) :: nomtyp(*)
    character(len=*) :: nofaex(*)
    character(len=80) :: nogrfa(nbgrou)
    character(len=*) :: nomamd
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE - FORMAT MED - LES FAMILLES - 2
!
! --------------------------------------------------------------------------------------------------
!
!     L'ENSEMBLE DES FAMILLES EST L'INTERSECTION DE L'ENSEMBLE
!     DES GROUPES : UN NOEUD/MAILLE APPARAIT AU PLUS DANS 1 FAMILLE
!     TABLE  NUMEROS DES FAMILLES POUR LES NOEUDS  <-> TABLE  DES COO
!     TABLES NUMEROS DES FAMILLES POUR MAILLE/TYPE <-> TABLES DES CNX
!     PAR CONVENTION, LES FAMILLES DE NOEUDS SONT NUMEROTEES >0 ET LES
!     FAMILLES DE MAILLES SONT NUMEROTEES <0. LA FAMILLE NULLE EST
!     DESTINEE AUX NOEUDS / ELEMENTS SANS FAMILLE.
! ENTREES :
!   FID    : IDENTIFIANT DU FICHIER MED
!   NOMAMD : NOM DU MAILLAGE MED
!   TYPENT : TYPE D'ENTITES : 0, POUR DES NOEUDS, 1 POUR DES MAILLES
!   NBRENT : NOMBRE D'ENTITES A TRAITER
!   NBGROU : NOMBRE DE GROUPES D'ENTITES
!   NOMGEN : VECTEUR NOMS DES GROUPES D'ENTITES
!   NBEC   : NOMBRE D'ENTIERS CODES
!   NOMAST : NOM DU MAILLAGE ASTER
!   PREFIX : PREFIXE POUR LES TABLEAUX DES RENUMEROTATIONS
!   TYPGEO : TYPE MED POUR CHAQUE MAILLE
!   NOMTYP : NOM DES TYPES POUR CHAQUE MAILLE
!   NMATYP : NOMBRE DE MAILLES PAR TYPE
! TABLEAUX DE TRAVAIL
!   NUFAEN : NUMERO DE FAMILLE POUR CHAQUE ENTITE
!            PAR DEFAUT, L'ALLOCATION AVEC JEVEUX A TOUT MIS A 0. CELA
!            SIGNIFIE QUE LES ENTITES APPARTIENNENT A LA FAMILLE NULLE.
!   NUFACR : NUMERO DE FAMILLES CREES. AU MAXIMUM, AUTANT QUE D'ENTITES
!   NOGRFA : NOM DES GROUPES ASSOCIES A CHAQUE FAMILLE.
!   NOFAEX = NOMS DES FAMILLES DEJA CREEES
!   TABAUX : PRESENCE D UNE ENTITE DANS UN GROUPE
! DIVERS
!   INFMED : NIVEAU DES INFORMATIONS SPECIFIQUES A MED A IMPRIMER
!   NIVINF : NIVEAU DES INFORMATIONS GENERALES
!   IFM    : UNITE LOGIQUE DU FICHIER DE MESSAGE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: edmail = 0, ednoeu = 3, tygeno = 0
    integer(kind=8) :: codret
    integer(kind=8) :: iaux, jaux, kaux
    integer(kind=8) :: numfam, nfam
    integer(kind=8) :: ityp
    integer(kind=8) :: nbeg, ige, ient, entfam, nbgnof, natt
    integer(kind=8) :: jgren
    integer(kind=8) :: tbaux(1)
    character(len=7) :: saux07
    character(len=8) :: saux08
    character(len=9) :: saux09
    character(len=64) :: nomfam
!
! --------------------------------------------------------------------------------------------------
!
!
!
!     NATT = NOMBRE D'ATTRIBUTS DANS UNE FAMILLE : JAMAIS. ELLES NE SONT
!            DEFINIES QUE PAR LES GROUPES
    natt = 0
!
!     NFAM = NUMERO DE LA DERNIERE FAMILLE ENREGISTREE (DE 0 A N>0)
!     FAMILLE 0 = ENTITES N'APPARTENANT A AUCUN GROUPE
    nfam = 0
!
!====
! 2. EN PRESENCE DE GROUPES, ON CREE DES FAMILLES
!====
!
    if (nbgrou .ne. 0) then
!
        if (typent .eq. tygeno) then
            saux09 = '.GROUPENO'
        else
            saux09 = '.GROUPEMA'
        end if
!
! 2.1. ==> BUT DE L'ETAPE 2.1 : CONNAITRE POUR CHAQUE ENTITE SES GROUPES
!          D'APPARTENANCE
!
        do ige = 1, nbgrou
            call jeveuo(jexnum(nomast//saux09, ige), 'L', jgren)
            call jelira(jexnum(nomast//saux09, ige), 'LONMAX', nbeg)
            if (infmed .ge. 2) then
                if (typent .eq. tygeno) then
                    saux07 = 'NOEUDS '
                else
                    saux07 = 'MAILLES'
                end if
                write (ifm, 100) nomgen(ige), nbeg, saux07
            end if
100         format('. GROUPE ', a, ' :', i12, 1x, a)
!           POUR CHAQUE GROUPE, ON BOUCLE SUR LES ENTITES QU'IL CONTIENT.
            do iaux = 1, nbeg
!
!           DEBUT VECTEUR ENTIER CODE POUR ENTITE IENT DANS JENTXG
                ient = zi(jgren-1+iaux)
                if (ient .ne. 0) then
!             ENREGISTREMENT APPARTENANCE DU ENTITE AU GROUPE
                    call setgfa(tabaux(1+(ient-1)*nbec), ige)
!             MISE A -1 DU NUM DE FAMILLE POUR CETTE ENTITE DANS NUFAEN
!             POUR INDIQUER QU'ELLE APPARTIENT AU MOINS A UN GROUPE
                    nufaen(ient) = 1
                end if
            end do
        end do
!
! 2.2. ==> BUT DE L'ETAPE 2.2 : FAIRE LA PARTITION EN FAMILLE ET NOTER :
!          . LE NUMERO DE LA 1ER ENTITE DE LA FAMILLE
!          . LE NUMERO DE FAMILLE DE CHAQUE ENTITE
!
!          ON BOUCLE SUR LES ENTITES APPARTENANT AU MOINS A UN GROUPE
!          ET ON LES RASSEMBLE PAR IDENTITE D'APPARTENANCE.
!          LES FAMILLES SONT NUMEROTEES DANS L'ORDRE D'APPARITION
!          ATTENTION : CET ALGORITHME A ETE OPTIMISE LE 6/9/2002
!                      ETRE PRUDENT DANS LES AMELIORATIONS FUTURES ...
!                      LES SITUATIONS PENALISANTES SONT CELLES-CI :
!                      QUELQUES DIZAINES DE MILLIERS D'ENTITES ET
!                      QUELQUES CENTAINES DE GROUPES
!                      EXEMPLE : ON EST PASSE D'UNE VINGTAINE D'HEURES
!                      A 3 MINUTES AVEC UN GROS MAILLAGE :
!                      . 426 817 NOEUDS EN 57 GROUPES ET
!                      . 418 514 MAILLES EN 8 629 GROUPES.
!
        do ient = 1, nbrent
            if (nufaen(ient) .ne. 0) then
!
!         BOUCLE 221 : ON PARCOURT TOUTES LES FAMILLES DEJA VUES.
!         POUR CHACUNE D'ELLES, ON COMPARE LES GROUPES ASSOCIES ET LES
!         GROUPES DE L'ENTITE COURANTE :
!         MEME COMPOSITION DE GROUPES <==> MEMES ENTIERS CODES
!         . SI C'EST LA MEME COMPOSITION DE GROUPES, LA FAMILLE EST LA
!           MEME. ON DONNE DONC LE NUMERO DE FAMILLE L'ENTITE COURANTE.
!         . SI ON N'A TROUVE AUCUNE FAMILLE, C'EST QU'UNE NOUVELLE
!           FAMILLE VIENT D'APPARAITRE. ON STOCKE SES CARACTERISTIQUES.
!
                jaux = nbec*(ient-1)
                do numfam = 1, nfam
                    entfam = nufacr(numfam)
                    kaux = nbec*(entfam-1)
                    do iaux = 1, nbec
                        if (tabaux(jaux+iaux) .ne. tabaux(kaux+iaux)) then
                            goto 221
                        end if
                    end do
!             ON A TROUVE UNE FAMILLE AVEC LA MEME COMPOSITION :
!             . ON NOTE QUE LA FAMILLE EST LA MEME
!             . ON PASSE A L'ENTITE SUIVANTE
                    nufaen(ient) = nufaen(entfam)
                    goto 22
221                 continue
                end do
!           AUCUN ENTITE NE CORRESPONDAIT : ON CREE UNE NOUVELLE FAMILLE
                nfam = nfam+1
!           ON MEMORISE CE NUMERO DE FAMILLE POUR L'ENTITE COURANTE
!           ATTENTION : LA CONVENTION MED VEUT QUE LE NUMERO SOIT
!           POSITIF POUR LES FAMILLES DE NOEUDS, NEGATIF POUR
!           LES MAILLES
                nufaen(ient) = nfam
                if (typent .ne. tygeno) then
                    nufaen(ient) = -nufaen(ient)
                end if
!           ON INDIQUE OU SE TROUVE LA 1ERE REFERENCE A CETTE FAMILLE
!           DANS LE VECTEUR NUFACR POUR EVITER DE PERDRE SON TEMPS APRES
                nufacr(nfam) = ient
            end if
22          continue
        end do
!
! 2.3. ==> BUT DE L'ETAPE 2.3 : CREATION DES FAMILLES D'ENTITES ET LES
!          ECRIRE DANS LE FICHIER
!
!          ON PARCOURT LES FAMILLES REPERTORIEES.
!          ON MEMORISE LES NOMS ET NUMEROS DES GROUPES QUI LA
!          CARACTERISENT. POUR CELA, ON SE BASE SUR LE PREMIER ENTITE
!          QUI EN FAIT PARTIE.
!
        do iaux = 1, nfam
!
! 2.3.1. ==> DETERMINATION DE LA FAMILLE : NOM, NOMS ET NUMEROS DES
!              GROUPES ASSOCIES
            numfam = iaux
            if (typent .ne. tygeno) then
                numfam = -numfam
            end if
!
!         NUMERO DE LA 1ERE ENTITE FAISANT REFERENCE A CETTE FAMILLE
            ient = nufacr(iaux)
!
!         NB ET NOMS+NUMS DES GROUPES ASSOCIES A LA FAMILLE
            call nomgfa(nomgen, nbgrou, tabaux(1+(ient-1)*nbec), nogrfa, nbgnof)
!
!         NOM DE LA FAMILLE : ON LE CONSTRUIT A PARTIR DES NOMS
!         DE GROUPES
!
            jaux = iaux-1
            call mdnofa(numfam, nogrfa, nbgnof, jaux, nofaex, nomfam)
!
! 2.3.2. ==> INFORMATION EVENTUELLE
!
            if (infmed .ge. 2) then
                jaux = 0
                do ient = 1, nbrent
                    if (nufaen(ient) .eq. numfam) then
                        jaux = jaux+1
                    end if
                end do
                if (typent .eq. tygeno) then
                    kaux = 0
                else
                    kaux = jaux
                    jaux = 0
                end if
                call desgfa(typent+1, numfam, nomfam, nbgnof, nogrfa, &
                            natt, tbaux, jaux, kaux, ifm, &
                            codret)
            end if
!
! 2.3.3. ==> ECRITURE DES CARACTERISTIQUES DE LA FAMILLE
!
            call as_mfacre(fid, nomamd, nomfam, numfam, nbgnof, nogrfa, codret)
            if (codret .ne. 0) then
                saux08 = 'mfacre'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        end do
    end if
!
!====
! 3. ECRITURE DE LA TABLE DES NUMEROS DE FAMILLES DES ENTITES
!    CELA SE FAIT PAR TYPE. ON REUTILISE LES VECTEURS CONTENANT
!    LES NUMEROS D'ENTITES/TYPE
!====
!
! 3.1. ==> ECRITURE DANS LE CAS DES NOEUDS
!
    if (typent .eq. tygeno) then
        call as_mmhfnw(fid, nomamd, nufaen, nbrent, ednoeu, tygeno, codret)
        if (codret .ne. 0) then
            saux08 = 'mmhfnw'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
! 3.2. ==> ECRITURE DANS LE CAS DES MAILLES : IL FAUT PASSER PAR LA
!          RENUMEROTATION ASTER-MED
!
    else
        do ityp = 1, MT_NTYMAX
            if (nmatyp(ityp) .ne. 0) then
!               RECUPERATION DU TABLEAU DES RENUMEROTATIONS
                call jeveuo('&&'//prefix//'.NUM.'//nomtyp(ityp), 'L', kaux)
!               CREATION VECTEUR NUMEROS DE FAMILLE POUR LES MAILLES / TYPE
                do iaux = 1, nmatyp(ityp)
                    tabaux(iaux) = nufaen(zi(kaux-1+iaux))
                end do
                call as_mmhfnw(fid, nomamd, tabaux, nmatyp(ityp), edmail, typgeo(ityp), codret)
                if (codret .ne. 0) then
                    saux08 = 'mmhfnw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
            end if
        end do
    end if
!
end subroutine
