! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine irmmf3(fid, nomamd, typent, nbrent, nbgrou,&
                  nomgen, nbec, nomast, prefix, typgeo,&
                  nomtyp, nmatyp, nufaen, nufacr, nogrfa,&
                  nofaex, tabaux, infmed, ifm, nosdfu)
!
implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_bcast_char80.h"
#include "asterc/asmpi_bcast_i.h"
#include "asterfort/as_mfacre.h"
#include "asterfort/as_mfrall.h"
#include "asterfort/as_mfrblc.h"
#include "asterfort/as_mfrdea.h"
#include "asterfort/as_mmhaaw.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/desgfa.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mdnofa.h"
#include "asterfort/nomgfa.h"
#include "asterfort/setgfa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
med_idt :: fid
integer :: typgeo(*), nmatyp(*)
integer :: typent, nbrent, nbgrou
integer :: nbec
integer :: nufaen(nbrent), nufacr(nbrent), tabaux(*)
integer :: infmed
integer :: ifm
character(len=6) :: prefix
character(len=8) :: nomast
character(len=24) :: nomgen(*)
character(len=8) :: nomtyp(*)
character(len=*) :: nofaex(*)
character(len=80) :: nogrfa(nbgrou)
character(len=*) :: nomamd
character(len=8) :: nosdfu
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
    character(len=6), parameter :: nompro = 'IRMMF3'
    integer, parameter :: edmail = 0, ednoeu = 3, tygeno = 0
    integer :: edfuin
    parameter (edfuin=0)
    integer :: codret
    integer :: iaux, jaux, kaux
    integer :: numfam, nfam, cmpt
    integer :: ityp, jnbno, jno, jma, nbnot, nbnol, start, filter(1)
    integer :: nbeg, ige, ient, entfam, nbgnof, natt, nbmal, nbmat, jtyp
    integer :: jgren, compt, jtest3, jtest4, jtest5, jtest6
    integer :: nbgr, nbgrp, nfam_max
    integer :: rang, nbproc, jgrou, jnufa, numgrp, jnofa, jnbgr, jtest, jtest2
    character(len=8) :: saux08
    character(len=9) :: saux09
    character(len=80) :: nomfam
    real(kind=8) :: start_time, end_time
    aster_logical :: lfamtr
    mpi_int :: mrank, msize, world, proc, taille
!
! --------------------------------------------------------------------------------------------------
!
!
    if (infmed .gt. 1) then
        call cpu_time(start_time)
        write (ifm,*) '<',nompro,'> DEBUT ECRITURE DES FAMILLES MED EN PARALLELE : '
    endif
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
    call asmpi_comm('GET', world)
    call asmpi_info(rank = mrank, size = msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    if (nbgrou .ne. 0) then
!
        if (typent .eq. tygeno) then
            saux09 = '.GROUPENO'
        else
            saux09 = '.GROUPEMA'
        endif
!
! 2.1. ==> BUT DE L'ETAPE 2.1 : CONNAITRE POUR CHAQUE ENTITE SES GROUPES
!          D'APPARTENANCE
!
        do ige = 1 , nbgrou
            call jeexin(jexnom(nomast//saux09, nomgen(ige)), codret)
            if(codret.ne.0) then
                call jeveuo(jexnom(nomast//saux09, nomgen(ige)), 'L', jgren)
                call jelira(jexnom(nomast//saux09, nomgen(ige)), 'LONMAX', nbeg)
!           POUR CHAQUE GROUPE, ON BOUCLE SUR LES ENTITES QU'IL CONTIENT.
                do iaux = 1 , nbeg
!
!           DEBUT VECTEUR ENTIER CODE POUR ENTITE IENT DANS JENTXG
                    ient = zi(jgren-1+iaux)
                    if (ient .ne. 0) then
!             ENREGISTREMENT APPARTENANCE DU ENTITE AU GROUPE
                        call setgfa(tabaux(1+(ient-1)*nbec), ige)
!             MISE A -1 DU NUM DE FAMILLE POUR CETTE ENTITE DANS NUFAEN
!             POUR INDIQUER QU'ELLE APPARTIENT AU MOINS A UN GROUPE
                        nufaen(ient) = 1
                    endif
                end do
            endif
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
                do numfam = 1 , nfam
                    entfam = nufacr(numfam)
                    kaux = nbec*(entfam-1)
                    do iaux = 1 , nbec
                        if (tabaux(jaux+iaux) .ne. tabaux(kaux+iaux)) then
                            goto 221
                        endif
                    end do
!             ON A TROUVE UNE FAMILLE AVEC LA MEME COMPOSITION :
!             . ON NOTE QUE LA FAMILLE EST LA MEME
!             . ON PASSE A L'ENTITE SUIVANTE
                    nufaen(ient) = nufaen(entfam)
                    goto 22
    221             continue
                end do
!           AUCUN ENTITE NE CORRESPONDAIT : ON CREE UNE NOUVELLE FAMILLE
                nfam = nfam + 1
!           ON MEMORISE CE NUMERO DE FAMILLE POUR L'ENTITE COURANTE
!           ATTENTION : LA CONVENTION MED VEUT QUE LE NUMERO SOIT
!           POSITIF POUR LES FAMILLES DE NOEUDS, NEGATIF POUR
!           LES MAILLES
                nufaen(ient) = nfam
                if (typent .ne. tygeno) then
                    nufaen(ient) = -nufaen(ient)
                endif
!           ON INDIQUE OU SE TROUVE LA 1ERE REFERENCE A CETTE FAMILLE
!           DANS LE VECTEUR NUFACR POUR EVITER DE PERDRE SON TEMPS APRES
                nufacr(nfam) = ient
            endif
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
        nfam_max = nfam
        call asmpi_comm_vect('MPI_MAX', 'I', sci=nfam_max)
        if(nfam_max.ne.0) then
            call wkvect('&&IRMMF2.NOGRFA', 'V V K80', nbgrou*nbgrou, jgrou)
            call wkvect('&&IRMMF2.NOMFAM', 'V V K80', nfam_max, jnofa)
            call wkvect('&&IRMMF2.NBGRFA', 'V V I', nfam_max, jnbgr)
            call wkvect('&&IRMMF2.NUMFAM', 'V V I', nfam_max, jnufa)
            numgrp = 0
            do iaux = 1 , nfam
!
! 2.3.1. ==> DETERMINATION DE LA FAMILLE : NOM, NOMS ET NUMEROS DES
!              GROUPES ASSOCIES
                numfam = iaux
                if (typent .ne. tygeno) then
                    numfam = -numfam
                endif
!
!         NUMERO DE LA 1ERE ENTITE FAISANT REFERENCE A CETTE FAMILLE
                ient = nufacr(iaux)
!
!         NB ET NOMS+NUMS DES GROUPES ASSOCIES A LA FAMILLE
                codret = tabaux(1+(ient-1)*nbec)
                call nomgfa(nomgen, nbgrou, tabaux(1+(ient-1)*nbec), nogrfa, nbgnof)
                zi(jnufa+iaux-1) = tabaux(1+(ient-1)*nbec)
!
!         NOM DE LA FAMILLE : ON LE CONSTRUIT A PARTIR DES NOMS
!         DE GROUPES
!
                jaux = iaux - 1
                call mdnofa(numfam, nogrfa, nbgnof, jaux, nofaex, nomfam)
                zi(jnbgr+iaux-1) = nbgnof
                zk80(jnofa+iaux-1) = nomfam
                do jaux = 1, nbgnof
                    zk80(jgrou+numgrp) = nogrfa(jaux)
                    numgrp = numgrp + 1
                enddo
            end do
            compt = 0
            call wkvect('&&IRMMF2.TEST', 'V V I', nbproc, jtest)
            zi(jtest+rang) = nfam
            call asmpi_comm_vect('MPI_SUM', 'I', nbval=nbproc, vi=zi(jtest))
            do jaux = 0, nbproc-1
                compt = compt + zi(jtest+jaux)
            enddo
            call wkvect('&&IRMMF2.TEST3', 'V V K80', compt, jtest3)
            call wkvect('&&IRMMF2.TEST4', 'V V I', compt, jtest4)
            compt = 0
            do jaux = 0, nbproc-1
                call wkvect('&&IRMMF2.TEST2', 'V V K80', max(1,zi(jtest+jaux)), jtest2)
                call wkvect('&&IRMMF2.TEST5', 'V V I', max(1,zi(jtest+jaux)), jtest5)
                if(jaux.eq.rang) then
                    do iaux = 1, zi(jtest+jaux)
                        zk80(jtest2+iaux-1) = zk80(jnofa+iaux-1)
                        zi(jtest5+iaux-1) = zi(jnbgr+iaux-1)
                    enddo
                endif
                proc = to_mpi_int(jaux)
                taille = to_mpi_int(zi(jtest+jaux))
                call asmpi_bcast_char80(zk80(jtest2), taille, proc, world)
                call asmpi_bcast_i(zi(jtest5), taille, proc, world)

                nbgr = 0
                do iaux = 1, taille
                    nbgr = nbgr + zi(jtest5+iaux-1)
                enddo
                call wkvect('&&IRMMF2.TEST6', 'V V K80', max(1,nbgr), jtest6)
                if(jaux.eq.rang) then
                    nbgrp = 0
                    do iaux = 1, nfam
                        do kaux = 1, zi(jnbgr+iaux-1)
                            zk80(jtest6+nbgrp) = zk80(jgrou+nbgrp)
                            nbgrp = nbgrp + 1
                        enddo
                    enddo
                    ASSERT(nbgrp.eq.nbgr)
                endif
                taille = to_mpi_int(nbgr)
                call asmpi_bcast_char80(zk80(jtest6), taille, proc, world)
!
                nbgr = 0
                do iaux = 1, zi(jtest+jaux)
                    lfamtr = .false._1
                    do kaux = 1, compt
                        if(zk80(jtest3+kaux-1).eq.zk80(jtest2+iaux-1)) then
                            lfamtr = .true._1
                            exit
                        endif
                    enddo
                    if(.not.lfamtr) then
                        zk80(jtest3+compt) = zk80(jtest2+iaux-1)
                        compt = compt + 1
                        nomfam = zk80(jtest2+iaux-1)
                        nbgnof = zi(jtest5+iaux-1)
                        do kaux = 1, nbgnof
                            nogrfa(kaux) = zk80(jtest6+nbgr+kaux-1)
                        enddo
!
! 2.3.2. ==> ECRITURE DES CARACTERISTIQUES DE LA FAMILLE
!
                        if (typent .ne. tygeno) then
                            call as_mfacre(fid, nomamd, nomfam, -compt, nbgnof, nogrfa, codret)
                        else
                            call as_mfacre(fid, nomamd, nomfam, compt, nbgnof, nogrfa, codret)
                        endif
                        if (codret .ne. 0) then
                            saux08='mfacre'
                            call utmess('F', 'DVP_97', sk=saux08, si=codret)
                        endif
                        if(jaux.eq.rang) then
                            zi(jtest4+iaux-1) = compt
                        endif
                    else
                        if(jaux.eq.rang) then
                            zi(jtest4+iaux-1) = kaux
                        endif
                    endif
                    nbgr = nbgr + zi(jtest5+iaux-1)
                enddo
                call jedetr('&&IRMMF2.TEST2')
                call jedetr('&&IRMMF2.TEST5')
                call jedetr('&&IRMMF2.TEST6')
            enddo
            call jedetr('&&IRMMF2.NOGRFA')
            call jedetr('&&IRMMF2.NOMFAM')
            call jedetr('&&IRMMF2.NBGRFA')
            call jedetr('&&IRMMF2.NUMFAM')
            call jedetr('&&IRMMF2.TEST')
            call jedetr('&&IRMMF2.TEST3')
        endif
        if (typent .ne. tygeno) then
            do ient = 1, nbrent
                if(nufaen(ient).eq.0) then
                    nufaen(ient) = 0
                else
                    nufaen(ient) = -zi(jtest4-nufaen(ient)-1)
                endif
            enddo
        else
            call jeveuo(nosdfu//'.NOEU', 'L', jno)
            cmpt = 0
            do ient = 1, nbrent
                if(zi(jno+ient-1).gt.0) then
                    cmpt = cmpt + 1
                    if(nufaen(ient).eq.0) then
                        nufaen(cmpt) = 0
                    else
                        nufaen(cmpt) = zi(jtest4+nufaen(ient)-1)
                    endif
                endif
            enddo
        endif
        call jedetr('&&IRMMF2.TEST4')
    endif
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
        call jeveuo(nosdfu//'.NBNO', 'L', jnbno)
        start = zi(jnbno)
        nbnol = zi(jnbno+1)
        nbnot = zi(jnbno+2)
        ASSERT(cmpt.eq.nbnol)
        call as_mfrall(1, filter, codret)
!
        call as_mfrblc(fid, nbnot, 1, 1, 0,&
                       edfuin, 2, "", start, nbnol,&
                       1, nbnol, 0, filter(1), codret)
        if (codret .ne. 0) then
            saux08='mfrblc'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        endif
        call as_mmhaaw(fid, nomamd, nufaen, nbnol, filter(1),&
                       ednoeu, tygeno, codret)
        if (codret .ne. 0) then
            saux08='mmhaaw'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        endif
!
        call as_mfrdea(1, filter, codret)
        if (codret .ne. 0) then
            saux08='mfrdea'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        endif
!
! 3.2. ==> ECRITURE DANS LE CAS DES MAILLES : IL FAUT PASSER PAR LA
!          RENUMEROTATION ASTER-MED
!
    else
        call jeveuo(nosdfu//'.MAIL', 'L', jma)
        call jeveuo(nosdfu//'.MATY', 'L', jtyp)
        do ityp = 1 , MT_NTYMAX
            if (zi(jtyp+3*(ityp-1)+2).ne.0) then
                start = zi(jtyp+3*(ityp-1))
                nbmal = nmatyp(ityp)
                nbmat = zi(jtyp+3*(ityp-1)+2)
                call as_mfrall(1, filter, codret)
                call as_mfrblc(fid, nbmat, 1, 1, 0,&
                               edfuin, 2, "", start, nbmal,&
                               1, nbmal, 0, filter(1), codret)
                if (codret .ne. 0) then
                    saux08='mfrblc'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                endif
!
!               RECUPERATION DU TABLEAU DES RENUMEROTATIONS
                call jeveuo('&&'//prefix//'.NUM.'//nomtyp(ityp), 'L', kaux)
!               CREATION VECTEUR NUMEROS DE FAMILLE POUR LES MAILLES / TYPE
                cmpt = 0
                do iaux = 1 , nmatyp(ityp)
                    if(zi(jma+ient-1).gt.0) then
                        cmpt = cmpt + 1
                        tabaux(cmpt) = nufaen(zi(kaux-1+iaux))
                    endif
                end do
                ASSERT(cmpt.eq.nbmal)
!
                call as_mmhaaw(fid, nomamd, tabaux, nbmal, filter(1),&
                            edmail, typgeo(ityp), codret)
                if (codret .ne. 0) then
                    saux08='mmhaaw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                endif
!
                call as_mfrdea(1, filter, codret)
                if (codret .ne. 0) then
                    saux08='mfrdea'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                endif
            endif
        enddo
    endif
!
    if (infmed .gt. 1) then
        call cpu_time(end_time)
        write (ifm,*) '<',nompro,'> FIN ECRITURE DES FAMILLES MED EN PARALLELE EN ', &
            end_time-start_time, "sec."
    endif
!
end subroutine
