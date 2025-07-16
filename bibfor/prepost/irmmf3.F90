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
subroutine irmmf3(fid, nomamd, typent, nbrent, nbgrou, &
                  nomgen, nbec, nomast, prefix, typgeo, &
                  nomtyp, nmatyp, nufaen, nufacr, nogrfa, &
                  nofaex, tabaux, infmed, ifm, nosdfu)
!
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_allgatherv_i.h"
#include "asterc/asmpi_allgatherv_char80.h"
#include "asterc/asmpi_allgather_i.h"
#include "asterfort/as_mfacre.h"
#include "asterfort/as_mfrall.h"
#include "asterfort/as_mfrblc.h"
#include "asterfort/as_mfrdea.h"
#include "asterfort/as_mmhaaw.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mdnofa.h"
#include "asterfort/nomgfa.h"
#include "asterfort/setgfa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
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
    character(len=8) :: nosdfu
!
! -----------------------------------------------------------------------------------------------
!
!     ECRITURE DU MAILLAGE - FORMAT MED - LES FAMILLES - 2
!
! -----------------------------------------------------------------------------------------------
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
! -----------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRMMF3'
    integer(kind=8), parameter :: edmail = 0, ednoeu = 3, tygeno = 0
    integer(kind=8) :: edfuin
    parameter(edfuin=0)
    integer(kind=8) :: codret
    integer(kind=8) :: iaux, jaux, kaux
    integer(kind=8) :: numfam, nfam, cmpt, ii
    integer(kind=8) :: ityp, jnbno, jno, jma, nbnot, nbnol, start, filter(1)
    integer(kind=8) :: nbeg, ige, ient, entfam, nbgnof, natt, nbmal, nbmat, jtyp
    integer(kind=8) :: jgren, jtest4, i_fama, kfama
    integer(kind=8) :: nbgr, nfam_max, nbbloc, nbfam_tot, nbgr_tot
    integer(kind=8) :: rang, nbproc, jgrou, jnufa, numgrp, jnofa, jnbgr, jtest, jtest12
    character(len=8) :: saux08
    character(len=9) :: saux09
    character(len=64) :: nomfam
    aster_logical :: lfamtr
    real(kind=8) :: start_time, end_time, start1, end1, start2, end2
    mpi_int :: mrank, msize, world, taille, one4
    character(len=80), pointer :: v_nomfag(:) => null()
    character(len=80), pointer :: v_nomgfag(:) => null()
    character(len=80), pointer :: v_fama(:) => null()
    integer(kind=8), pointer :: v_nbgrg(:) => null()
    mpi_int, pointer :: v_count(:) => null()
    mpi_int, pointer :: v_displ(:) => null()
    mpi_int, pointer :: v_count2(:) => null()
    mpi_int, pointer :: v_displ2(:) => null()

!
! ----------------------------------------------------------------------------------------------
!
!
    if (typent .eq. tygeno) then
        saux09 = '.GROUPENO'
        saux08 = "NOEUDS"
    else
        saux09 = '.GROUPEMA'
        saux08 = "MAILLES"
    end if

    if (infmed .gt. 1) then
        call cpu_time(start_time)
        write (ifm, *) '<', nompro, '> DEBUT ECRITURE DES FAMILLES DE ' &
            //saux08//' MED EN PARALLELE : '
    end if
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
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    if (nbgrou .ne. 0) then
!
!
! 2.1. ==> BUT DE L'ETAPE 2.1 : CONNAITRE POUR CHAQUE ENTITE SES GROUPES
!          D'APPARTENANCE
!
        call cpu_time(start1)
        do ige = 1, nbgrou
            call jeexin(jexnom(nomast//saux09, nomgen(ige)), codret)
            if (codret .ne. 0) then
                call jeveuo(jexnom(nomast//saux09, nomgen(ige)), 'L', jgren)
                call jelira(jexnom(nomast//saux09, nomgen(ige)), 'LONMAX', nbeg)
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
            end if
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
        call cpu_time(end1)
        if (infmed .gt. 1) then
            write (ifm, *) '<', nompro, '> ** Preparation des familles en ', end1-start1, "sec."
        end if
!
! 2.3. ==> BUT DE L'ETAPE 2.3 : CREATION DES FAMILLES D'ENTITES ET LES
!          ECRIRE DANS LE FICHIER
!
!          ON PARCOURT LES FAMILLES REPERTORIEES.
!          ON MEMORISE LES NOMS ET NUMEROS DES GROUPES QUI LA
!          CARACTERISENT. POUR CELA, ON SE BASE SUR LE PREMIER ENTITE
!          QUI EN FAIT PARTIE.
!
        call cpu_time(start1)
        nfam_max = nfam
        call asmpi_comm_vect('MPI_MAX', 'I', sci=nfam_max)
        if (nfam_max .ne. 0) then
            call cpu_time(start2)
!
            ! On cherche le nombre de groupe en local
            numgrp = 0
            ! boucle sur famille locale
            do iaux = 1, nfam
!
!         NUMERO DE LA 1ERE ENTITE FAISANT REFERENCE A CETTE FAMILLE
                ient = nufacr(iaux)
!
!         NB ET NOMS+NUMS DES GROUPES ASSOCIES A LA FAMILLE
                call nomgfa(nomgen, nbgrou, tabaux(1+(ient-1)*nbec), nogrfa, nbgnof)
                numgrp = numgrp+nbgnof
            end do
!
            call wkvect('&&IRMMF2.NOGRFA', 'V V K80', max(1, numgrp), jgrou)
            call wkvect('&&IRMMF2.NOMFAM', 'V V K80', max(1, nfam), jnofa)
            call wkvect('&&IRMMF2.NBGRFA', 'V V I', max(1, nfam), jnbgr)
            call wkvect('&&IRMMF2.NUMFAM', 'V V I', max(1, nfam), jnufa)
!
            call wkvect('&&IRMMF2.TEST', 'V V I', nbproc, jtest)
            call wkvect('&&IRMMF2.COUNT', 'V V S', nbproc, vi4=v_count)
            call wkvect('&&IRMMF2.DISPL', 'V V S', nbproc+1, vi4=v_displ)
!
            one4 = to_mpi_int(1)
            call asmpi_allgather_i([nfam], one4, zi(jtest), one4, world)
!
            nbfam_tot = 0
            ! on compte le nb de famille totale
            do jaux = 0, nbproc-1
                nbfam_tot = nbfam_tot+zi(jtest+jaux)
                v_count(jaux+1) = to_mpi_int(zi(jtest+jaux))
                v_displ(jaux+2) = to_mpi_int(nbfam_tot)
            end do
!
! On a plus de famille en HPC qu'en std car les familles sont redécoupés par sous-domaine
! donc 1 famille en std peut donner jusqu'à nb_sous_domaine familles dans le maillage de fin
!
            if (infmed .gt. 1) then
                write (ifm, *) '<', nompro, '> ** Nombre de familles locales/globales: ', &
                    nfam, nbfam_tot, numgrp
            end if
!
            numfam = 0
            do jaux = 0, rang
                numfam = numfam+zi(jtest+jaux)
            end do
!
            numgrp = 0
            ! boucle sur famille locale
            do iaux = 1, nfam
!
! 2.3.1. ==> DETERMINATION DE LA FAMILLE : NOM, NOMS ET NUMEROS DES
!              GROUPES ASSOCIES
                numfam = numfam+iaux
                if (typent .ne. tygeno) then
                    numfam = -numfam
                end if
!
!         NUMERO DE LA 1ERE ENTITE FAISANT REFERENCE A CETTE FAMILLE
                ient = nufacr(iaux)
!
!         NB ET NOMS+NUMS DES GROUPES ASSOCIES A LA FAMILLE
                codret = tabaux(1+(ient-1)*nbec)
                call nomgfa(nomgen, nbgrou, tabaux(1+(ient-1)*nbec), nogrfa, nbgnof)
                zi(jnufa+iaux-1) = tabaux(1+(ient-1)*nbec)
!
!         NOM DE LA FAMILLE : ON LE CONSTRUIT A PARTIR DU NUMERO DE FAMILLE
!
                jaux = iaux-1
                call mdnofa(numfam, nogrfa, nbgnof, jaux, nofaex, nomfam)
!
                ! nb de groupe de la famille
                zi(jnbgr+iaux-1) = nbgnof
                ! nom de la famille
                zk80(jnofa+iaux-1) = nomfam
                ! noms des groupes de la famille
                do jaux = 1, nbgnof
                    zk80(jgrou+numgrp) = nogrfa(jaux)
                    numgrp = numgrp+1
                end do
            end do
!
            call cpu_time(end2)
            if (infmed .gt. 1) then
                write (ifm, *) '<', nompro, '> ** Création des ' &
                    //'noms de familles en ', end2-start2, ' sec'
            end if
            call cpu_time(start1)
!
            call wkvect('&&IRMMF2.TEST4', 'V V I', max(1, nfam), jtest4)
            call wkvect('&&IRMMF2.NBGRG', 'V V I', nbfam_tot, vi=v_nbgrg)
            call wkvect('&&IRMMF2.NOMFAG', 'V V K80', nbfam_tot, vk80=v_nomfag)
!
            ! AllGather des noms et nb groupe par famille
            taille = to_mpi_int(nfam)
            call asmpi_allgatherv_i(zi(jnbgr), taille, v_nbgrg, v_count, v_displ, world)
            call asmpi_allgatherv_char80(zk80(jnofa), taille, v_nomfag, v_count, v_displ, world)

            call wkvect('&&IRMMF2.TEST12', 'V V I', nbproc, jtest12)
            call wkvect('&&IRMMF2.COUNT2', 'V V S', nbproc, vi4=v_count2)
            call wkvect('&&IRMMF2.DISPL2', 'V V S', nbproc+1, vi4=v_displ2)
            one4 = to_mpi_int(1)
            call asmpi_allgather_i([numgrp], one4, zi(jtest12), one4, world)
!
            ! on compte le nb de grp
            nbgr_tot = 0
            do jaux = 0, nbproc-1
                nbgr_tot = nbgr_tot+zi(jtest12+jaux)
                v_count2(jaux+1) = to_mpi_int(zi(jtest12+jaux))
                v_displ2(jaux+2) = to_mpi_int(nbgr_tot)
            end do

            call wkvect('&&IRMMF2.NOMGFAG', 'V V K80', nbgr_tot, vk80=v_nomgfag)
            ! Allgather des groupes
            taille = to_mpi_int(numgrp)
            call asmpi_allgatherv_char80(zk80(jgrou), taille, v_nomgfag, v_count2, &
                                         v_displ2, world)
!
            call cpu_time(end1)
            if (infmed .gt. 1) then
                write (ifm, *) '<', nompro, '> ** Partage des familles en ', &
                    end1-start1, "sec."
            end if
!
            ! boucle sur les familles
            call cpu_time(start1)
            call wkvect('&&IRMMF2.FAMA', 'V V K80', nbfam_tot, vk80=v_fama)
            nbgr = 0
            i_fama = 0
            ii = 0
            do iaux = 1, nbfam_tot
                nomfam = v_nomfag(iaux)
                nbgnof = v_nbgrg(iaux)
!
                lfamtr = ASTER_FALSE
                do kaux = 1, i_fama
                    if (v_fama(kaux) .eq. nomfam) then
                        lfamtr = ASTER_TRUE
                        kfama = kaux
                        exit
                    end if
                end do
!
                if (.not. lfamtr) then
                    do kaux = 1, nbgnof
                        nogrfa(kaux) = v_nomgfag(nbgr+kaux)
                    end do
!
! 2.3.2. ==> ECRITURE DES CARACTERISTIQUES DE LA FAMILLE
!
                    i_fama = i_fama+1
                    if (typent .ne. tygeno) then
                        call as_mfacre(fid, nomamd, nomfam, -i_fama, nbgnof, nogrfa, codret)
                    else
                        call as_mfacre(fid, nomamd, nomfam, i_fama, nbgnof, nogrfa, codret)
                    end if
                    if (codret .ne. 0) then
                        saux08 = 'mfacre'
                        call utmess('F', 'DVP_97', sk=saux08, si=codret)
                    end if
                    kfama = i_fama
                    v_fama(i_fama) = nomfam
                end if
                nbgr = nbgr+nbgnof
!
                if (iaux > v_displ(rang+1) .and. iaux <= v_displ(rang+2)) then
                    ii = ii+1
                    zi(jtest4+ii-1) = kfama
                end if
            end do
            ASSERT(ii == nfam)
!
            call jedetr('&&IRMMF2.NOGRFA')
            call jedetr('&&IRMMF2.NOMFAM')
            call jedetr('&&IRMMF2.NBGRFA')
            call jedetr('&&IRMMF2.NUMFAM')
            call jedetr('&&IRMMF2.TEST')
            call jedetr('&&IRMMF2.TEST12')
            call jedetr('&&IRMMF2.COUNT')
            call jedetr('&&IRMMF2.DISPL')
            call jedetr('&&IRMMF2.COUNT2')
            call jedetr('&&IRMMF2.DISPL2')
            call jedetr('&&IRMMF2.NBGRG')
            call jedetr('&&IRMMF2.NOMFAG')
            call jedetr('&&IRMMF2.NOMGFAG')
            call jedetr('&&IRMMF2.FAMA')
        end if
        if (typent .ne. tygeno) then
            do ient = 1, nbrent
                if (nufaen(ient) .eq. 0) then
                    nufaen(ient) = 0
                else
                    nufaen(ient) = -zi(jtest4-nufaen(ient)-1)
                end if
            end do
        else
            call jeveuo(nosdfu//'.NOEU', 'L', jno)
            cmpt = 0
            do ient = 1, nbrent
                if (zi(jno+ient-1) .gt. 0) then
                    cmpt = cmpt+1
                    if (nufaen(ient) .eq. 0) then
                        nufaen(cmpt) = 0
                    else
                        nufaen(cmpt) = zi(jtest4+nufaen(ient)-1)
                    end if
                end if
            end do
        end if
        call jedetr('&&IRMMF2.TEST4')
    end if
    call cpu_time(end1)
    if (infmed .gt. 1) then
        write (ifm, *) '<', nompro, '> ** Ecriture caractéristiques des familles en ', &
            end1-start1, "sec."
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
    call cpu_time(start1)
    if (typent .eq. tygeno) then
        call jeveuo(nosdfu//'.NBNO', 'L', jnbno)
        start = zi(jnbno)
        nbnol = zi(jnbno+1)
        nbnot = zi(jnbno+2)
        ASSERT(cmpt .eq. nbnol)
        call as_mfrall(1, filter, codret)
!
        call as_mfrblc(fid, nbnot, 1, 1, 0, &
                       edfuin, 2, "", start, nbnol, &
                       1, nbnol, 0, filter(1), codret)
        if (codret .ne. 0) then
            saux08 = 'mfrblc'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
        call as_mmhaaw(fid, nomamd, nufaen, nbnol, filter(1), &
                       ednoeu, tygeno, codret)
        if (codret .ne. 0) then
            saux08 = 'mmhaaw'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
        call as_mfrdea(1, filter, codret)
        if (codret .ne. 0) then
            saux08 = 'mfrdea'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
! 3.2. ==> ECRITURE DANS LE CAS DES MAILLES : IL FAUT PASSER PAR LA
!          RENUMEROTATION ASTER-MED
!
    else
        call jeveuo(nosdfu//'.MAIL', 'L', jma)
        call jeveuo(nosdfu//'.MATY', 'L', jtyp)
        do ityp = 1, MT_NTYMAX
            nbmat = zi(jtyp+3*(ityp-1)+2)
            if (nbmat .ne. 0) then
                start = zi(jtyp+3*(ityp-1))
                nbmal = nmatyp(ityp)
                ASSERT(zi(jtyp+3*(ityp-1)+1) == nbmal)
                nbbloc = 1
                if (nbmal == 0) nbbloc = 0
                call as_mfrall(1, filter, codret)
                call as_mfrblc(fid, nbmat, 1, 1, 0, &
                               edfuin, 2, "", start, nbmal, &
                               nbbloc, nbmal, 0, filter(1), codret)
                if (codret .ne. 0) then
                    saux08 = 'mfrblc'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                cmpt = 0
                if (nbmal > 0) then
!               RECUPERATION DU TABLEAU DES RENUMEROTATIONS
                    call jeveuo('&&'//prefix//'.NUM.'//nomtyp(ityp), 'L', kaux)
!               CREATION VECTEUR NUMEROS DE FAMILLE POUR LES MAILLES / TYPE
                    do iaux = 1, nmatyp(ityp)
                        if (zi(jma+ient-1) .gt. 0) then
                            cmpt = cmpt+1
                            tabaux(cmpt) = nufaen(zi(kaux-1+iaux))
                        end if
                    end do
                end if
                ASSERT(cmpt .eq. nbmal)
!
                call as_mmhaaw(fid, nomamd, tabaux, nbmal, filter(1), &
                               edmail, typgeo(ityp), codret)
                if (codret .ne. 0) then
                    saux08 = 'mmhaaw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                call as_mfrdea(1, filter, codret)
                if (codret .ne. 0) then
                    saux08 = 'mfrdea'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
            end if
        end do
    end if
    call cpu_time(end1)
    if (infmed .gt. 1) then
        write (ifm, *) '<', nompro, '> ** Ecriture des familles en ', end1-start1, "sec."
    end if
!
    if (infmed .gt. 1) then
        call cpu_time(end_time)
        write (ifm, *) '<', nompro, '> FIN ECRITURE DES FAMILLES DE ' &
            //saux08//' MED EN PARALLELE EN ', end_time-start_time, "sec."
    end if
!
end subroutine
