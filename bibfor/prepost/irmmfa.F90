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
subroutine irmmfa(fid, nomamd, nbnoeu, nbmail, nomast, &
                  nbgrno, nomgno, nbgrma, nomgma, prefix, &
                  typgeo, nomtyp, nmatyp, infmed, nosdfu)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/as_mfacre.h"
#include "asterfort/infniv.h"
#include "asterfort/irmmf1.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    med_idt :: fid
    integer(kind=8) :: typgeo(*), nmatyp(*)
    integer(kind=8) :: nbnoeu, nbmail, nbgrno, nbgrma
    integer(kind=8) :: infmed
    character(len=6) :: prefix
    character(len=8) :: nomast
    character(len=24) :: nomgno(*), nomgma(*)
    character(len=8) :: nomtyp(*)
    character(len=*) :: nomamd
    character(len=8) :: nosdfu
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE - FORMAT MED - LES FAMILLES
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
!     ENTREE:
!       FID    : IDENTIFIANT DU FICHIER MED
!       NOMAMD : NOM DU MAILLAGE MED
!       NBNOEU : NOMBRE DE NOEUDS DU MAILLAGE
!       NBMAIL : NOMBRE DE MAILLES DU MAILLAGE
!       NOMAST : NOM DU MAILLAGE ASTER
!       NBGRNO : NOMBRE DE GROUPES DE NOEUDS
!       NBGRMA : NOMBRE DE GROUPES DE MAILLES
!       NOMGNO : VECTEUR NOMS DES GROUPES DE NOEUDS
!       NOMGMA : VECTEUR NOMS DES GROUPES DE MAILLES
!       PREFIX : PREFIXE POUR LES TABLEAUX DES RENUMEROTATIONS
!       TYPGEO : TYPE MED POUR CHAQUE MAILLE
!       NOMTYP : NOM DES TYPES POUR CHAQUE MAILLE
!       NMATYP : NOMBRE D'ENTITES PAR TYPE
!       INFMED : NIVEAU DES INFORMATIONS A IMPRIMER
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRMMFA'
    integer(kind=8) :: tygeno
    integer(kind=8) :: codret
    integer(kind=8) :: iaux
    integer(kind=8) :: numfam
    integer(kind=8) :: natt
    integer(kind=8) :: jnofam
    integer(kind=8) :: jmafam
    integer(kind=8) :: ifm, niv
    character(len=8) :: saux08
    character(len=24) :: nufano, nufama
    character(len=64) :: nomfam
    character(len=80) :: saux80
    real(kind=8) :: start_time, end_time
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call infniv(ifm, niv)
    if (niv .gt. 1) then
        call cpu_time(start_time)
        write (ifm, *) '<', nompro, '> DEBUT ECRITURE DES FAMILLES MED : '
    end if
!
!     VECTEUR NUMEROS DES FAMILLES DES ENTITES = NB ENTITES
!     PAR DEFAUT, JEVEUX MET TOUT A 0. CELA SIGNIFIE QUE, PAR DEFAUT,
!     LES ENTITES APPARTIENNENT A LA FAMILLE NULLE.
!               12   345678   9012345678901234
    tygeno = 0
    nufano = '&&'//nompro//'.NU_FAM_NOE     '
    nufama = '&&'//nompro//'.NU_FAM_MAI     '
!
    call wkvect(nufano, 'V V I', nbnoeu, jnofam)
    call wkvect(nufama, 'V V I', nbmail, jmafam)
!
!====
! 2. LES FAMILLES DE NOEUDS
!====
!
!
    iaux = tygeno
    call irmmf1(fid, nomamd, iaux, nbnoeu, nbgrno, &
                nomgno, zi(jnofam), nomast, prefix, typgeo, &
                nomtyp, nmatyp, infmed, ifm, nosdfu)
!
!GN      WRITE (IFM,*)
!GN     >'TEMPS CPU POUR CREER/ECRIRE LES FAMILLES DE NOEUDS  :',TPS2(4)
!
!====
! 3. LES FAMILLES DE MAILLES
!====
!
!
    iaux = tygeno+1
    call irmmf1(fid, nomamd, iaux, nbmail, nbgrma, &
                nomgma, zi(jmafam), nomast, prefix, typgeo, &
                nomtyp, nmatyp, infmed, ifm, nosdfu)
!
!GN      WRITE (IFM,*)
!GN     >'TEMPS CPU POUR CREER/ECRIRE LES FAMILLES DE MAILLES :',TPS2(4)
!
!====
! 4. ON CREE TOUJOURS UNE FAMILLE DE NUMERO 0 NE REFERENCANT RIEN,
!    POUR LES NOEUDS ET ELEMENTS N'APPARTENANT A AUCUN GROUPE
!    REMARQUE : IL FAUT LE FAIRE A LA FIN POUR AVOIR LES BONNES
!    IMRPESSIONS
!====
!
! 4.1. ==> CARACTERISTIQUE
!
    numfam = 0
    natt = 0
    nomfam = 'FAMILLE_NULLE___________________'//'________________________________'
!
! 4.2. ==> ECRITURE
!
    call as_mfacre(fid, nomamd, nomfam, numfam, 0, &
                   saux80, codret)
    if (codret .ne. 0) then
        saux08 = 'mfacre'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
!====
! 5. LA FIN
!====
!
    call jedetr(nufano)
    call jedetr(nufama)
!
!GN      WRITE (IFM,*) '==> DUREE TOTALE DE ',NOMPRO,' :',TPS1(4)
    if (niv .gt. 1) then
        call cpu_time(end_time)
        write (ifm, *) '<', nompro, '> FIN ECRITURE DES FAMILLES MED EN ', &
            end_time-start_time, "sec."
    end if
!
    call jedema()
!
end subroutine
