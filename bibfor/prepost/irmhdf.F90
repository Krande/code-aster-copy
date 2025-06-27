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
subroutine irmhdf(ifi, ndim, nbnoeu, coordo, nbmail, &
                  connex, point, nomast, typma, titre, &
                  nbtitr, nbgrno, nomgno, nbgrma, nomgma, &
                  nommai, nomnoe, infmed, lfichUniq, nosdfu)
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mmhcre.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/irmdes.h"
#include "asterfort/irmmfa.h"
#include "asterfort/irmmma.h"
#include "asterfort/irmmno.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jemarq.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/mdexma.h"
#include "asterfort/mdnoma.h"
#include "asterfort/ulisog.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/irmhpc.h"
!
    integer(kind=8) :: connex(*), typma(*), point(*)
    integer(kind=8) :: ifi, ndim, nbnoeu, nbmail, nbgrno, nbgrma
    integer(kind=8) :: infmed, nbtitr
    character(len=80) :: titre(*)
    character(len=8) :: nommai(*), nomnoe(*), nomast
    character(len=24) :: nomgno(*), nomgma(*)
    real(kind=8) :: coordo(*)
    aster_logical, optional :: lfichUniq
    character(len=8), optional :: nosdfu
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE - FORMAT MED
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREE:
!       IFI    : UNITE LOGIQUE D'IMPRESSION DU MAILLAGE
!       NDIM   : DIMENSION DU PROBLEME (2  OU 3)
!       NBNOEU : NOMBRE DE NOEUDS DU MAILLAGE
!       COORDO : VECTEUR DES COORDONNEES DES NOEUDS
!       NBMAIL : NOMBRE DE MAILLES DU MAILLAGE
!       CONNEX : CONNECTIVITES
!       POINT  : VECTEUR POINTEUR DES CONNECTIVITES (LONGUEURS CUMULEES)
!       NOMAST : NOM DU MAILLAGE
!       TYPMA  : VECTEUR TYPES DES MAILLES
!       TITRE  : TITRE ASSOCIE AU MAILLAGE
!       NBGRNO : NOMBRE DE GROUPES DE NOEUDS
!       NBGRMA : NOMBRE DE GROUPES DE MAILLES
!       NOMGNO : VECTEUR NOMS DES GROUPES DE NOEUDS
!       NOMGMA : VECTEUR NOMS DES GROUPES DE MAILLES
!       NOMMAI : VECTEUR NOMS DES MAILLES
!       NOMNOE : VECTEUR NOMS DES NOEUDS
!       INFMED : NIVEAU DES INFORMATIONS A IMPRIMER
!       LFICHUNIQ: ASTER LOGICAL, FICHIER UNIQUE
!       NOSDFU   : NOM STRUCTURE DONNEE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRMHDF'
    integer(kind=8), parameter :: edlect = 0, edleaj = 1, edcrea = 3, ednstr = 0, edcart = 0
    integer(kind=8) :: edmode, codret
    integer(kind=8) :: nbtyp
    med_idt :: fid, ifimed
    integer(kind=8) :: nmatyp(MT_NTYMAX), nnotyp(MT_NTYMAX), typgeo(MT_NTYMAX)
    integer(kind=8) :: renumd(MT_NTYMAX), modnum(MT_NTYMAX), numnoa(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: iaux, jaux, nuanom(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: lnomam
    integer(kind=8) :: ifm, niv
    real(kind=8) :: start_time, end_time
    character(len=1) :: saux01
    character(len=6) :: saux06
    character(len=8) :: nomtyp(MT_NTYMAX)
    character(len=8) :: saux08, nom_sd_par
    character(len=16) :: saux16(0:3)
    character(len=64) :: nomamd
    character(len=80) :: descdt
    character(len=200) :: nofimd, desc
    character(len=255) :: kfic
    character(len=64) :: valk(2)
    aster_logical :: existm, ficexi, lfu
    character(len=16), parameter :: nocoor(3) = (/'X               ', &
                                                  'Y               ', &
                                                  'Z               '/)
    character(len=16), parameter :: uncoor(3) = (/'INCONNU         ', &
                                                  'INCONNU         ', &
                                                  'INCONNU         '/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call infniv(ifm, niv)
    if (niv .gt. 1) then
        call cpu_time(start_time)
        write (ifm, *) '<', nompro, '> DEBUT ECRITURE DU MAILLAGE MED : '
    end if
!
    if (present(lfichUniq)) then
        nom_sd_par = ' '
        if (lfichUniq) then
            ASSERT(present(nosdfu))
            nom_sd_par = nosdfu
        end if
        lfu = lfichUniq
    else
        nom_sd_par = ' '
        lfu = .false._1
    end if
!
! 1.2. ==> NOM DU FICHIER MED
!
    call ulisog(ifi, kfic, saux01)
    if (kfic(1:1) .eq. ' ') then
        call codent(ifi, 'G', saux08)
        nofimd = 'fort.'//saux08
    else
        nofimd = kfic(1:200)
    end if
    if (lfu .and. nofimd(1:1) .ne. "/") then
        call utmess("F", "MED_11", sk=nofimd)
    end if
!
    if (niv .gt. 1) then
        write (ifm, *) '<', nompro, '> NOM DU FICHIER MED : ', nofimd
    end if
!
! - Generate name of mesh for MED
!
    call mdnoma(nomamd, lnomam, nomast, codret)
!
!
! 1.4. ==> LE MAILLAGE EST-IL DEJA PRESENT DANS LE FICHIER ?
!          SI OUI, ON NE FAIT RIEN DE PLUS QU'EMETTRE UNE INFORMATION
!
    iaux = 0
    ifimed = 0
    call mdexma(nofimd, ifimed, nomamd, iaux, existm, &
                jaux, codret)
!
    if (existm) then
!
        valk(1) = nofimd(1:32)
        valk(2) = nomamd
        call utmess('A', 'MED_67', nk=2, valk=valk)
!
!     ------------------------------------------------------------------
!
    else
!
!====
! 2. DEMARRAGE
!====
!
! 2.1. ==> OUVERTURE FICHIER MED EN MODE
!      SOIT 'CREATION' SI LE FICHIER N'EXISTE PAS ENCORE,
!      SOIT 'LECTURE_AJOUT' (CELA SIGNIFIE QUE LE FICHIER EST ENRICHI).
!
!     TEST L'EXISTENCE DU FICHIER
        inquire (file=nofimd, exist=ficexi)
        if (ficexi) then
            edmode = edlect
            if (lfu) then
                call as_med_open(fid, nofimd, edmode, codret, .true._1)
            else
                call as_med_open(fid, nofimd, edmode, codret)
            end if
            if (codret .ne. 0) then
                edmode = edcrea
            else
                edmode = edleaj
                call as_mficlo(fid, codret)
                if (codret .ne. 0) then
                    saux08 = 'mficlo'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
            end if
        else
            edmode = edcrea
        end if
        if (lfu) then
            call as_med_open(fid, nofimd, edmode, codret, .true._1)
        else
            call as_med_open(fid, nofimd, edmode, codret)
        end if
        if (codret .ne. 0) then
            saux08 = 'mfiope'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
        if (infmed .ge. 2) then
            saux16(edlect) = 'LECTURE SEULE.  '
            saux16(edleaj) = 'LECTURE/ECRITURE'
            saux16(edcrea) = 'CREATION.       '
            call codent(edmode, 'G', saux08)
            valk(1) = saux08
            valk(2) = saux16(edmode)
            call utmess('I', 'MED_40', nk=2, valk=valk)
        end if
!
! 2.2. ==> CREATION DU MAILLAGE AU SENS MED (TYPE MED_NON_STRUCTURE)
!
!GN      PRINT *,'APPEL DE as_mmhcre AVEC :'
!GN      PRINT *,NOMAMD
!GN      PRINT *,NDIM
!GN      PRINT *,EDNSTR
        desc = 'CREE PAR CODE_ASTER'
        descdt = 'SANS UNITES'
        call as_mmhcre(fid, nomamd, ndim, ednstr, desc, &
                       descdt, edcart, nocoor, uncoor, codret)
        if (codret .ne. 0) then
            saux08 = 'mmhcre'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
! 2.3. ==> . RECUPERATION DES NB/NOMS/NBNO/NBITEM DES TYPES DE MAILLES
!            DANS CATALOGUE
!          . RECUPERATION DES TYPES GEOMETRIE CORRESPONDANT POUR MED
!          . VERIF COHERENCE AVEC LE CATALOGUE
!
        call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, &
                    modnum, nuanom, numnoa)
!
!====
! 3. LA DESCRIPTION
!====
!
        if (edmode .eq. edcrea) then
!
            call irmdes(fid, titre, nbtitr, infmed)
!
        end if
!
!====
! 4. LES NOEUDS
!====
!
        call irmmno(fid, nomamd, ndim, nbnoeu, coordo, &
                    nomnoe, nom_sd_par)
!
!====
! 5. LES MAILLES
!====
!
        saux06 = nompro
!
        call irmmma(fid, nomamd, nbmail, connex, point, &
                    typma, nommai, saux06, nbtyp, typgeo, &
                    nomtyp, nnotyp, renumd, nmatyp, infmed, &
                    modnum, nuanom, nom_sd_par)
!
!====
! 6. LES FAMILLES
!====
!
        saux06 = nompro
!
        call irmmfa(fid, nomamd, nbnoeu, nbmail, nomast, &
                    nbgrno, nomgno, nbgrma, nomgma, saux06, &
                    typgeo, nomtyp, nmatyp, infmed, nom_sd_par)
!
!====
! 7. LES EQUIVALENCES
!====
!
!     CALL IRMMEQ ()  ! NE FAIT RIEN ...
!
!====
! 8. IMPRESSION NUMEROTATION GLOBALE ET JOINTS EN HPC
!====
!
        if (isParallelMesh(nomast) .and. (.not. lfu)) then
            call irmhpc(fid, nomamd, nomast, nbnoeu)
        end if
!
!====
! 9. FERMETURE DU FICHIER MED
!====
!
        call as_mficlo(fid, codret)
        if (codret .ne. 0) then
            saux08 = 'mficlo'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
!====
! 10. LA FIN
!====
!
        call jedetc('V', '&&'//nompro, 1)
!
    end if
!
    if (niv .gt. 1) then
        call cpu_time(end_time)
        write (ifm, *) '<', nompro, '> FIN ECRITURE DU MAILLAGE MED EN ', &
            end_time-start_time, "sec."
    end if
!
    call jedema()
!
end subroutine
