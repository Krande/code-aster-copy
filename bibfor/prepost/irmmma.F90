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
subroutine irmmma(fid, nomamd, nbCell, connex, point, &
                  typma, nommai, prefix, nbtyp, typgeo, &
                  nomtyp, nnotyp, renumd, nbCellType, infmed, &
                  modnum, nuanom, nosdfu)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/as_mfrall.h"
#include "asterfort/as_mfrblc.h"
#include "asterfort/as_mfrdea.h"
#include "asterfort/as_mmhcyw.h"
#include "asterfort/as_mmhenw.h"
#include "asterfort/as_mmhyaw.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    med_idt :: fid
    integer(kind=8) :: nbCell, nbtyp
    integer(kind=8) :: connex(*), typma(*), point(*)
    integer(kind=8) :: typgeo(*), nnotyp(*), nbCellType(MT_NTYMAX)
    integer(kind=8) :: renumd(*), modnum(MT_NTYMAX), nuanom(MT_NTYMAX, *)
    integer(kind=8) :: infmed
    character(len=6) :: prefix
    character(len=8) :: nommai(*)
    character(len=8) :: nomtyp(*)
    character(len=*) :: nomamd
    character(len=8) :: nosdfu
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE -  FORMAT MED - LES MAILLES
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREE:
!       FID    : IDENTIFIANT DU FICHIER MED
!       NOMAMD : NOM DU MAILLAGE MED
!       NBMAIL : NOMBRE DE MAILLES DU MAILLAGE
!       CONNEX : CONNECTIVITES
!       POINT  : VECTEUR POINTEUR DES CONNECTIVITES (LONGUEURS CUMULEES)
!       TYPMA  : VECTEUR TYPES DES MAILLES
!       NOMMAI : VECTEUR NOMS DES MAILLES
!       PREFIX : PREFIXE POUR LES TABLEAUX DES RENUMEROTATIONS
!                A UTILISER PLUS TARD
!       NBTYP  : NOMBRE DE TYPES POSSIBLES POUR MED
!       TYPGEO : TYPE MED POUR CHAQUE MAILLE
!       NNOTYP : NOMBRE DE NOEUDS POUR CHAQUE TYPE DE MAILLES
!       NOMTYP : NOM DES TYPES POUR CHAQUE MAILLE
!       RENUMD : RENUMEROTATION DES TYPES ENTRE MED ET ASTER
!       INFMED : NIVEAU DES INFORMATIONS A IMPRIMER
!       MODNUM : INDICATEUR SI LA SPECIFICATION DE NUMEROTATION DES
!                NOEUDS DES MAILLES EST DIFFERENTES ENTRE ASTER ET MED:
!                     MODNUM = 0 : NUMEROTATION IDENTIQUE
!                     MODNUM = 1 : NUMEROTATION DIFFERENTE
!       NUANOM : TABLEAU DE CORRESPONDANCE DES NOEUDS MED/ASTER.
!                NUANOM(ITYP,J): NUMERO DANS ASTER DU J IEME NOEUD DE LA
!                MAILLE DE TYPE ITYP DANS MED.
!       NOSDFU : NOM STRUCTURE DONNÉE ?
!
!     SORTIE:
!       NMATYP : NOMBRE DE MAILLES PAR TYPE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRMMMA'
    integer(kind=8) :: edfuin, edmail, ednoda
    integer(kind=8) :: codret
    integer(kind=8) :: ipoin, iCellType, letype
    integer(kind=8) :: ino, iCell, nbCellTypeTotal, nbbloc
    integer(kind=8) :: jnomma(MT_NTYMAX), jnumma(MT_NTYMAX), jcnxma(MT_NTYMAX)
    integer(kind=8) :: ifm, niv, jma, rang, nbproc, jtyp, nbmat, nbmal, start, jno
    integer(kind=8) :: filter(1), numno
    mpi_int :: mrank, msize
    character(len=8) :: saux08
    aster_logical :: lnocen, lfu
    real(kind=8) :: start1, end1
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
    if (niv .gt. 1) then
        call cpu_time(start1)
        write (ifm, *) '<', nompro, '> DEBUT ECRITURE DES MAILLES : '
    end if
!
!====
! 2. PREPARATION DES TABLEAUX PAR TYPE DE MAILLE
!====
!
! 2.1. ==> DECOMPTE DU NOMBRE DE MAILLES PAR TYPE
!          EN FAIT, ON VEUT JUSTE SAVOIR S'IL Y EN A OU PAS.
!
    lfu = .false._1
    if (nosdfu .ne. ' ') then
        lfu = .true._1
        call jeveuo(nosdfu//'.MAIL', 'L', jma)
        call jeveuo(nosdfu//'.NOEU', 'L', jno)
        call jeveuo(nosdfu//'.MATY', 'L', jtyp)
    else
        jma = 0
        jno = 0
        jtyp = 0
    end if
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    edfuin = 0
    edmail = 0
    ednoda = 0
    nbCellType(1:) = 0
    if (lfu) then
        do iCell = 1, nbCell
            if (zi(jma+iCell-1) .ne. 0) nbCellType(typma(iCell)) = nbCellType(typma(iCell))+1
        end do
    else
        do iCell = 1, nbCell
            nbCellType(typma(iCell)) = nbCellType(typma(iCell))+1
        end do
    end if
!
!     ON TRAITE LE TETRA15 -> TETRA10, ON OUBLIE LES 5 DERNIERS NOEUDS
!     ON L'IMPRIME COMME UN TETRA10
!     PYRAM19 -> PYRAM13, PENTA21 -> PENTA18, HEXA9 -> HEXA8 et PENTA7 -> PENTA6
    lnocen = ASTER_FALSE
    if (nbCellType(MT_TETRA15) .ne. 0) then
        nbCellType(MT_TETRA10) = nbCellType(MT_TETRA10)+nbCellType(MT_TETRA15)
        nbCellType(MT_TETRA15) = 0
        lnocen = ASTER_TRUE
    end if
    if (nbCellType(MT_PYRAM19) .ne. 0) then
        nbCellType(MT_PYRAM13) = nbCellType(MT_PYRAM13)+nbCellType(MT_PYRAM19)
        nbCellType(MT_PYRAM19) = 0
        lnocen = ASTER_TRUE
    end if
    if (nbCellType(MT_PENTA21) .ne. 0) then
        nbCellType(MT_PENTA18) = nbCellType(MT_PENTA18)+nbCellType(MT_PENTA21)
        nbCellType(MT_PENTA21) = 0
        lnocen = ASTER_TRUE
    end if
    if (nbCellType(MT_HEXA9) .ne. 0) then
        nbCellType(MT_HEXA8) = nbCellType(MT_HEXA8)+nbCellType(MT_HEXA9)
        nbCellType(MT_HEXA9) = 0
        lnocen = ASTER_TRUE
    end if
    if (nbCellType(MT_PENTA7) .ne. 0) then
        nbCellType(MT_PENTA6) = nbCellType(MT_PENTA6)+nbCellType(MT_PENTA7)
        nbCellType(MT_PENTA7) = 0
        lnocen = ASTER_TRUE
    end if
!     LES NOEUDS N'ONT PAS ÉTE RETIRÉS DU MAILLAGE, ILS SONT IMPRIMÉS AVEC IRMMNO
    if (lnocen) then
        call utmess('I', 'PREPOST_86')
    end if
!
! 2.2. ==> ON VERIFIE QUE L'ON SAIT ECRIRE LES MAILLES PRESENTES DANS
!          LE MAILLAGE
!
    do iCellType = 1, MT_NTYMAX
        if (nbCellType(iCellType) .ne. 0) then
            if (typgeo(iCellType) .eq. 0) then
                call utmess('F', 'PREPOST2_93', sk=nomtyp(iCellType))
            end if
        end if
    end do
!
! 2.3. ==> CREATION DE PLUSIEURS VECTEURS PAR TYPE DE MAILLE PRESENT :
!              UN VECTEUR CONTENANT LES NOMS DES MAILLES/TYPE
!           +  UN VECTEUR CONTENANT LES NUMEROS DES MAILLES/TYPE
!           +  UN VECTEUR CONTENANT LA CONNECTIVITE DES MAILLE/TYPE
!              (CONNECTIVITE = NOEUDS + UNE VALEUR BIDON(0) SI BESOIN)
!
    do iCellType = 1, MT_NTYMAX
        if (nbCellType(iCellType) .ne. 0) then
            call wkvect('&&'//nompro//'.NOM.'//nomtyp(iCellType), 'V V K16', &
                        nbCellType(iCellType), jnomma(iCellType))
            call wkvect('&&'//prefix//'.NUM.'//nomtyp(iCellType), 'V V I', &
                        nbCellType(iCellType), jnumma(iCellType))
            call wkvect('&&'//nompro//'.CNX.'//nomtyp(iCellType), 'V V I', &
                        nnotyp(iCellType)*nbCellType(iCellType), jcnxma(iCellType))
            ! nnotyp(iCellType) ASSOCIE MOINS DE NOEUDS AUX ELEMENTS NON SUPPORTES PAR MED
        end if
    end do
!
! 2.4. ==> ON PARCOURT TOUTES LES MAILLES. POUR CHACUNE D'ELLES, ON
!          STOCKE SON NOM, SON NUMERO, SA CONNECTIVITE
!          LA CONNECTIVITE EST FOURNIE EN STOCKANT TOUS LES NOEUDS A
!          LA SUITE POUR UNE MAILLE DONNEE.
!          C'EST CE QU'ON APPELLE LE MODE ENTRELACE DANS MED
!          A LA FIN DE CETTE PHASE, NMATYP CONTIENT LE NOMBRE DE MAILLES
!          POUR CHAQUE TYPE
!
    nbCellType = 0
!
    do iCell = 1, nbCell
        if (lfu) then
            if (zi(jma+iCell-1) .eq. 0) cycle
        end if
        iCellType = typma(iCell)
!       CAS PARTICULIER - MAILLE NON SUPPORTEE PAR MED
        if (iCellType .eq. MT_HEXA9) iCellType = MT_HEXA8
        if (iCellType .eq. MT_TETRA15) iCellType = MT_TETRA10
        if (iCellType .eq. MT_PYRAM19) iCellType = MT_PYRAM13
        if (iCellType .eq. MT_PENTA21) iCellType = MT_PENTA18
        if (iCellType .eq. MT_PENTA7) iCellType = MT_PENTA6
        ipoin = point(iCell)
        nbCellType(iCellType) = nbCellType(iCellType)+1
!       NOM DE LA MAILLE DE TYPE ITYP DANS VECT NOM MAILLES
        zk16(jnomma(iCellType)-1+nbCellType(iCellType)) = nommai(iCell)
!       NUMERO ASTER DE LA MAILLE DE TYPE ITYP DANS VECT NUM MAILLES
        zi(jnumma(iCellType)-1+nbCellType(iCellType)) = iCell
!       CONNECTIVITE DE LA MAILLE TYPE ITYP DANS VECT CNX:
!       I) POUR LES TYPES DE MAILLE DONT LA NUMEROTATION DES NOEUDS
!          ENTRE ASTER ET MED EST IDENTIQUE:
        if (modnum(iCellType) .eq. 0) then
            if (lfu) then
                do ino = 1, nnotyp(iCellType)
                    numno = abs(zi(jno+connex(ipoin-1+ino)-1))
                    zi(jcnxma(iCellType)-1+(nbCellType(iCellType)-1)*nnotyp(iCellType)+ino) = numno
                end do
            else
                do ino = 1, nnotyp(iCellType)
                    zi(jcnxma(iCellType)-1+(nbCellType(iCellType)-1)*nnotyp(iCellType)+ino) = &
                        connex(ipoin-1+ino)
                end do
            end if
!       II) POUR LES TYPES DE MAILLE DONT LA NUMEROTATION DES NOEUDS
!          ENTRE ASTER ET MED EST DIFFERENTE (CF LRMTYP):
        else
            if (lfu) then
                do ino = 1, nnotyp(iCellType)
                    numno = abs(zi(jno+connex(ipoin-1+nuanom(iCellType, ino))-1))
                    zi(jcnxma(iCellType)-1+(nbCellType(iCellType)-1)*nnotyp(iCellType)+ino) = numno
                end do
            else
                do ino = 1, nnotyp(iCellType)
                    zi(jcnxma(iCellType)-1+(nbCellType(iCellType)-1)*nnotyp(iCellType)+ino) = &
                        connex(ipoin-1+nuanom(iCellType, ino))
                end do
            end if
        end if
    end do
!
!====
! 3. ECRITURE
!    ON PARCOURT TOUS LES TYPES POSSIBLES POUR MED ET ON DECLENCHE LES
!    ECRITURES SI DES MAILLES DE CE TYPE SONT PRESENTES DANS LE MAILLAGE
!    LA RENUMEROTATION PERMET D'ECRIRE LES MAILLES DANS L'ORDRE
!    CROISSANT DE LEUR TYPE MED. CE N'EST PAS OBLIGATOIRE CAR ICI ON
!    FOURNIT LES TABLEAUX DE NUMEROTATION DES MAILLES. MAIS QUAND CES
!    TABLEAUX SONT ABSENTS, C'EST LA LOGIQUE QUI PREVAUT. DONC ON LA
!    GARDE DANS LA MESURE OU CE N'EST PAS PLUS CHER ET QUE C'EST CE QUI
!    EST FAIT A LA LECTURE
!====
!
    do letype = 1, nbtyp
!
! 3.0. ==> PASSAGE DU NUMERO DE TYPE MED AU NUMERO DE TYPE ASTER
!
        iCellType = renumd(letype)
!
        if (infmed .ge. 2) then
            write (ifm, 300) nomtyp(iCellType), nbCellType(iCellType)
        end if
300     format('TYPE ', a8, ' : ', i10, ' MAILLES')
!
        if (lfu) then
            nbCellTypeTotal = zi(jtyp+3*(iCellType-1)+2)
        else
            nbCellTypeTotal = nbCellType(iCellType)
        end if
!
        if (nbCellTypeTotal .ne. 0) then
!
! 3.1. ==> LES CONNECTIVITES
!          LA CONNECTIVITE EST FOURNIE EN STOCKANT TOUS LES NOEUDS A
!          LA SUITE POUR UNE MAILLE DONNEE.
!          C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!
            if (lfu) then
                nbmal = nbCellType(iCellType)
                ASSERT(nbmal .eq. zi(jtyp+3*(iCellType-1)+1))
                nbmat = zi(jtyp+3*(iCellType-1)+2)
                start = zi(jtyp+3*(iCellType-1))
                nbbloc = 1
                if (nbmal .eq. 0) nbbloc = 0
                call as_mfrall(1, filter, codret)
                call as_mfrblc(fid, nbmat, 1, nnotyp(iCellType), 0, &
                               0, 2, "", start, nbmal, &
                               nbbloc, nbmal, 0, filter(1), codret)
                if (codret .ne. 0) then
                    saux08 = 'mfrblc'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                call as_mmhyaw(fid, nomamd, zi(jcnxma(iCellType)), &
                               nnotyp(iCellType)*nbCellType(iCellType), edmail, &
                               typgeo(iCellType), ednoda, filter(1), codret)
                if (codret .ne. 0) then
                    saux08 = 'mmhcyw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                call as_mfrdea(1, filter, codret)
                if (codret .ne. 0) then
                    saux08 = 'mfrdea'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
            else
!
! 3.2. ==> LE NOM DES MAILLES
!
                call as_mmhcyw(fid, nomamd, zi(jcnxma(iCellType)), &
                               nnotyp(iCellType)*nbCellType(iCellType), edfuin, &
                               nbCellType(iCellType), edmail, typgeo(iCellType), ednoda, codret)
                if (codret .ne. 0) then
                    saux08 = 'mmhcyw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
! 3.3. ==> LE NUMERO DES MAILLES
!
                call as_mmhenw(fid, nomamd, zi(jnumma(iCellType)), nbCellType(iCellType), edmail, &
                               typgeo(iCellType), codret)
                if (codret .ne. 0) then
                    saux08 = 'mmhenw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
            end if
        end if
!
    end do
!
!====
! 4. LA FIN
!====
!
    do iCellType = 1, MT_NTYMAX
        if (nbCellType(iCellType) .ne. 0) then
            call jedetr('&&'//nompro//'.NOM.'//nomtyp(iCellType))
            call jedetr('&&'//nompro//'.CNX.'//nomtyp(iCellType))
        end if
    end do
!
    call jedema()
!
    if (niv .gt. 1) then
        call cpu_time(end1)
        write (ifm, *) '<', nompro, '> FIN ECRITURE DES MAILLES EN ', end1-start1, 'SEC'
    end if
!
end subroutine
