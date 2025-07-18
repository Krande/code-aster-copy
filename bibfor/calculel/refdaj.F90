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

subroutine refdaj(arret, result, nbordr, numer, typre, &
                  conre, codret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/refdag.h"
#include "asterfort/wkvect.h"
! person_in_charge: hassan.berro at edf.fr
!     ------------------------------------------------------------------
!                              FONCTION
!     _______________________________________________________________
!    | AJOUTER UNE NOUVELLE ENTREE DE REFERENCE A UN RESULTAT        |
!    | DYNAMIQUE : CALCUL VIBRATOIRE OU CALCUL MODAL                 |
!    |_______________________________________________________________|
!
! ---------
! EXEMPLES: call refdaj ('F',modmec,nbordr,numeq,'DYNAMIQUE'  ,matric,iret)
! --------- call refdaj ('F',tragen,nbordr,numddl,'DYNAMIQUE'  ,matric,iret)
!           call refdaj ('F',ritzba,nbordr,numddl,'INTERF_DYNA',intdyn,iret)
!           call refdaj ('F',ritzba,nbordr,numddl,'INTERF_STAT',modsta,iret)
!
!                     DESCRIPTIVE DES VARIABLES
!   ___________________________________________________________________
!  | IN > ARRET  : 'F' OU ' ' : ON ARRETE AVEC ERREUR FATALE EN    [K1]|
!  |               CAS DE PROBLEME DE COMPATIBILITE ?                  |
!  |                                                                   |
!  | IN > RESULT : RESULTAT DYNAMIQUE A ENRICHIR                   [K8]|
!  |OUT <                                                              |
!  |                                                                   |
!  | IN > NBORDR : NOMBRE DE NUMEROS D'ORDRE (CHAMPS) QUI SERONT   [I] |
!  |               ARCHIVES DANS LA SD DYNAMIQUE AVEC CES REFERENCES   |
!  |                                                                   |
!  | IN > NUMER  : NUMEROTATION DES CONCEPTS DE REFERENCE          [K*]|
!  |               SOIT : NUME_DDL                                     |
!  |               SOIT : NUME_EQUA                                    |
!  |                                                                   |
!  | IN > TYPRE  : TYPE DE REFERENCE A AJOUTER                     [K4]|
!  |               POSSIBILITES : 'DYNAMIQUE'                          |
!  |                              'INTERF_DYNA'                        |
!  |                              'INTERF_STAT'                        |
!  |                              'MESURE'                             |
!  |                                                                   |
!  | IN > CONREF : LISTE DES NOMS DE SD A REFERENCER          ARR [K24]|
!   ___________________________________________________________________
!  |OUT < CODRET : CODE RETOUR                                      [I]|
!  |               SOIT : > 0 SI L'OPERATION S'EST DEROULEE AVEC SUCCES|
!  |                      = 1 SI UNE MODIFICATION A ETE APPORTEE       |
!  |                      = 2 SI AUCUNE MODIFICATION A ETE NECESSAIRE  |
!  |                      = 0 SI L'OPERATION A ECHOUEE                 |
!   ___________________________________________________________________
!
!   ___________________________________________________________________
!
!  - 0 - INITIALISATIONS DIVERSES
!   ___________________________________________________________________
!
!     0.1 - DECLARATION DES VARIABLES D'ENTREE/SORTIE
!
    character(len=1) :: arret
    character(len=8) :: result
    character(len=*) :: numer
    character(len=*) :: typre
    character(len=*) :: conre(3)
    integer(kind=8) :: nbordr, codret
!
!     0.2 - DECLARATION DES VARIABLES LOCALES
!
    aster_logical :: oktres, newref, oktref
    integer(kind=8) :: lonref(4), indref, jrefe, nbrefs, nbrefsmax, nbinit, nbord1
    integer(kind=8) :: ibid, jbid, jindi, nbord0, ir, nbcham
    character(len=1) :: jvb
    character(len=8) :: k8bid, resu2
    character(len=24) :: typres, accres(10), accref(5), obindi, corefd, typref, kbid
    character(len=24) :: numer1, bl24, conref(3)
!
    data accres/'ACOU_HARMO', 'DYNA_HARMO', 'DYNA_TRANS', 'HARM_GENE', 'MODE_ACOU', &
        'MODE_FLAMB', 'MODE_MECA', 'MODE_MECA_C', 'TRAN_GENE', 'EVOL_NOLI'/
!
    data accref/'DYNAMIQUE', 'INTERF_DYNA', 'INTERF_STAT', 'MESURE', 'INIT'/
    data lonref/3, 1, 1, 1/
!
    typref = typre
    conref = conre
    numer1 = numer
    codret = 0
    bl24 = '                        '
    jvb = 'G'
    typres = ' '
    nbinit = 4
    nbord1 = nbordr
!
!     0.3 - ACTUALISATION DE LA VALEUR DE LA MARQUE COURANTE
!
    call jemarq()
!  ____________________________________________________________________
!
!  - 1 - VERIFICATION DES PARAMETRES D'ENTREE
!  ____________________________________________________________________
!
!     1.1 - CONDITION D'ARRET
    if ((arret .ne. 'F') .and. (arret .ne. ' ')) then
        ASSERT(.false.)
    end if
!
    call getres(resu2, typres, kbid)
    if (result(1:2) .eq. '&&') then
        jvb = 'V'
    else if (result .eq. ' ') then
        result = resu2
    else if (result .ne. resu2) then
        call gettco(result, typres)
    end if
!
!     1.2 - PRESENCE DU .INDI (OBJET INT) ET .REFD (COLLECTION D'OBJ K8)
    newref = .false.
    obindi = result//'           .INDI'
    corefd = result//'           .REFD'
    call jeexin(obindi, ibid)
    call jeexin(corefd, jbid)
    if ((ibid*jbid) .eq. 0 .and. (ibid+jbid) .ne. 0) then
        if (arret .eq. 'F') then
            ASSERT(.false.)
        end if
        codret = 0
        goto 27
    end if
    if ((ibid+jbid) .eq. 0) newref = .true.
!
!     1.3 - TYPE DE RESULTAT A TRAITER, CAS D'UN CONCEPT RE-ENTRANT
    if (.not. (newref)) then
        oktres = .false.
        do ibid = 1, 10
            if (typres .eq. accres(ibid)) oktres = .true.
        end do
        if (.not. (oktres)) then
            if (arret .eq. 'F') then
                ASSERT(.false.)
            end if
            codret = 0
            goto 27
        end if
    end if
!
!     1.4 - TYPE DE REFERENCE A AJOUTER
    oktref = .false.
    do ibid = 1, 5
        if (typref .eq. accref(ibid)) then
            indref = ibid
            oktref = .true.
            goto 16
        end if
    end do
16  continue
    if (.not. (oktref)) then
        if (arret .eq. 'F') then
            ASSERT(.false.)
        end if
        codret = 0
        goto 27
    end if
!  ____________________________________________________________________
!
!  - 2 - AJOUT DES ENTREES DE REFERENCE DANS LES .INDI ET .REFD
!  ____________________________________________________________________
!     2.1 - CREATION DU .INDI ET .REFD SI BESOIN
    if (newref) then
        call wkvect(obindi, jvb//' V I', nbinit, jbid)
        call jecrec(corefd, jvb//' V K24', 'NU', 'CONTIG', 'CONSTANT', &
                    nbinit)
        call jeecra(corefd, 'LONT', nbinit*5, k8bid)
!       INITIALISATION DES ELEMENTS DE .INDI A (-100)
        do ibid = 1, nbinit
            zi(jbid+ibid-1) = -100
        end do
        if (typref .eq. 'INIT') then
            codret = 1
            goto 27
        end if
    end if
!
!     2.2 - RAJOUT SI BESOIN DE L'ENTREE DE REF. A LA COLLECTION REFD
    call jelira(corefd, 'NUTIOC', nbrefs, k8bid)
    call jelira(corefd, 'LONT', nbrefsmax, k8bid)
    nbrefsmax = nbrefsmax/5
!
!     2.2.1 - VERIFIER QUE LE NOMBRE MAX D'OBJETS DE LA COLLECTION REFD
!             N'A PAS DEJA ETE ETEINT, DOUBLER SA TAILLE LE CAS ECHEANT
    if (nbrefs .ge. nbrefsmax) then
        call refdag(result)
    end if
!
!     2.2.2 - VERIFIER SI L'ENTREE DE REFERENCE EXISTE DEJA
!             CRITERES DE LA VERIFICATION :
!             1. TYPE D'ENTREE IDENTIQUE
!             2. NUMEROTATION IDENTIQUE
!             3. PREMIERE ENTREE IDENTIQUE
    if (nbrefs .ge. 1) then
        do ibid = 1, nbrefs
            call jeveuo(jexnum(corefd, ibid), 'E', jrefe)
            if ((typre .eq. zk24(jrefe)) .and. (numer .eq. zk24(jrefe+1)) .and. &
                (conre(1) .eq. zk24(jrefe+2))) then
!               --- ALORS METTRE A JOUR L'ENTREE DU .REFD
                do jbid = 1, lonref(indref)
                    zk24(jrefe+jbid+1) = conref(jbid)
                end do
!               --- VERIFIER EGALEMENT ET METTRE A JOUR LE .INDI SI LE NOMBRE DE
!                 - NUMEROS D'ORDRES INITIAL N'A PAS ETE CORRECTEMENT RENSEIGNE
!                 - (N2 - N1 = -1 / APPEL A REFDAJ AVEC NBCHAM = -1 )
                if (ibid .eq. nbrefs) then
                    nbord0 = 0
                    call jeveuo(obindi, 'E', jindi)
                    if (ibid .gt. 1) nbord0 = zi(jindi+nbrefs-1)
                    if ((zi(jindi+nbrefs)-nbord0) .eq. -1) then
                        zi(jindi+nbrefs) = nbord0+nbord1
                    else
                        zi(jindi+nbrefs-1) = zi(jindi+nbrefs-1)+nbord1
                    end if
                end if

                codret = 1
                goto 27
            end if
        end do
    end if
!
!   2.2.2 - RAJOUTER UN OBJET A LA COLLECTION .REFD AVEC LES INFOS
!           DE CONREF ET M.A.J L'OBJET .INDI
!
!   --- CUMULER L'INDICE CORRESPONDANT A LA NOUVELLE REFERENCE
    nbord0 = 0
    call jeveuo(obindi, 'E', jindi)
    if (nbrefs .ge. 1) then
        nbord0 = zi(jindi+nbrefs-1)
    else
!       --- TRAITEMENT SPECIAL SI IL S'AGIT DE LA PREMIERE ENTREE INDI :
!           SI nbord1 EST EGALE A -1 ET UN NOMBRE NON NULL DES CHAMPS EST DEJA STOCKE
!           DANS LA SD : ON Y MET CE NOMBRE
        call dismoi('NB_CHAMPS', result, 'RESU_DYNA', repi=nbcham, arret='C', &
                    ier=ir)
        if ((ir .eq. 0) .and. (nbcham .ne. 0)) then
            nbord1 = nbcham
        end if
    end if
!
    zi(jindi+nbrefs) = nbord0+nbord1
!
!   --- RAJOUTER UNE ENTREE A LA COLLECTION .REFD
    call jecroc(jexnum(corefd, nbrefs+1))
    call jeveuo(jexnum(corefd, nbrefs+1), 'E', jrefe)
    zk24(jrefe) = typre
    zk24(jrefe+1) = numer1
    do ibid = 1, lonref(indref)
        zk24(jrefe+ibid+1) = conref(ibid)
    end do
    if (lonref(indref) .lt. 3) then
        do ibid = lonref(indref)+1, 3
            zk24(jrefe+ibid+1) = bl24
        end do
    end if
    codret = 1
!
27  continue
!
!   For debugging purposes only
!   call utimsd(6, 1, .false._1, .true._1, corefd,1, 'G')
!   call utimsd(6, 1, .false._1, .true._1, obindi,1, 'G')
!
    call jedema()
!
end subroutine
