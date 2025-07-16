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

subroutine cgnoor(mafour, nomail, motfac, iocc, nbmc, &
                  motcle, typmcl, typlig, nbma, ndorig, &
                  ndextr, typm, vecori)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterc/getexm.h"
#include "asterfort/i2extf.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: iocc, nbmc, nbma
    character(len=*) :: typlig
    character(len=24) :: mafour
    character(len=8) :: nomail, ndorig, ndextr, typm
    character(len=16) :: motcle(*), typmcl(*)
    character(len=*) :: motfac
    real(kind=8) :: vecori(3)
! person_in_charge: jacques.pellet at edf.fr
!-----------------------------------------------------------------------
!  CETTE ROUTINE EST UTILISEE DANS DEFI_GROUP ET DEFI_FOND_FISS
!
! FONCTION REALISEE:
!  * CREER LA LISTE DES NUMEROS DES SEGX FORMANT UNE LIGNE (MAFOUR,NBMA)
!  * CALCULER LE NOM DES 2 NOEUDS EXTREMITES DE MAFOUR: NDORIG ET NDEXTR
!    (NDORIG EST TOUJOURS CALCULE. MAIS PARFOIS NDEXTR=' ')
!  * CALCULER LE TYPE DES SEGMENTS (SEG2/SEG3/SEG4)
!  * CALCULER LE VECTEUR D'ORIENTATION SI NECESSAIRE (VECORI)
!
! CETTE FONCTION REALISE AUSSI CERTAINES VERIFICATIONS :
!  * TOUTES LES MAILLES DE MAFOUR SONT DES SEGX
!  * TOUTES LES MAILLES DE MAFOUR SONT DE MEME TYPE SEG2/3/4
!  * MAFOUR DEFINIT UNE LIGNE CONTINUE SANS BIFURCATION
!
! ATTENTION : LA LISTE MAFOUR N'EST PAS TRIEE
!
!     ENTREES:
!        MAFOUR     : NOM DE L'OBJET QUI CONTIENDRA LES MAILLES SEG
!        NOMAIL     : NOM DU MAILLAGE
!        MOTFAC     : MOT-CLE FACTEUR
!        IOCC       : OCCURENCE COURANTE DE MOTFAC
!        NBMC       : NOMBRE DE MOT-CLE SIMPLE
!        MOTCLE     : MOT-CLE SIMPLE TYPE MA OU GROUP_MA A TRAITER
!        TYPMCL     : TYPE D'ENTITE ENTREE SOUS LE MOT-CLE
!        TYPLIG     : TYPE DE LIGNE : 'FERME' (OU NON)
!                     SI TYPLIG.EQ.' ' ET SI NDORIG=NDEXTR (/=' '),
!                     ON CONSIDERE QUE TYPLIG='FERME'
!     SORTIES:
!        MAFOUR     : L'OBJET EST ALLOUE ET REMPLI
!        NBMA       : NOMBRE DE MAILLES CONSIDEREES
!        NDORIG     : NOM DU NOEUD ORIGINE
!        NDEXTR     : NOM DU NOEUD EXTREMITE
!        TYPM       : TYPE DE MAILLE DE LA LISTE MAFOUR (SEG2/SEG3/SEG4)
!        VECORI     : 3 COORDONNEES EVENTUELLEMENT DONNEES PAR VECT_ORIE
!                     (0,0,0) SINON
!-----------------------------------------------------------------------
!
!
    integer(kind=8) :: jmail, jtypm, iatyma
    integer(kind=8) :: ier, im, n1, n2, n3, nid, nig, nbnot
    integer(kind=8) :: nunori, trouv, ibid, in, nd
    integer(kind=8) :: iret, ima, nbmato
    integer(kind=8) :: jcour2
    character(len=8) :: k8b, nomma, typmp
    character(len=16) :: k16bid, nomcmd, orig
    character(len=24) :: conec, typp, mesmai, valk(2), nogrp
    aster_logical :: bug, erreur, lexist
    integer(kind=8), pointer :: noeud_apparies(:) => null()
    integer(kind=8), pointer :: noeuds_extrem(:) => null()
    integer(kind=8), pointer :: compteur(:) => null()
    integer(kind=8), pointer :: type_noeud(:) => null()
! DEB-------------------------------------------------------------------
    call jemarq()
!
    call getres(k8b, k16bid, nomcmd)
!
!     ------------------------------------------------------------------
!     INITIALISATION DE VARIABLES
!     ------------------------------------------------------------------
    mesmai = '&&CGNOOR.MES_MAILLES'
    conec = nomail//'.CONNEX         '
    typp = nomail//'.TYPMAIL        '
    call dismoi('NB_NO_MAILLA', nomail, 'MAILLAGE', repi=nbnot)
    call jeveuo(typp, 'L', iatyma)
!
!
!     ------------------------------------------------------------------
!     RECUPERATION DES MAILLES SOUS LES MOT-CLES "MOTCLE"
!       => MESMAI : NOMS DES MAILLES
!       => MAFOUR : NUMEROS DES MAILLES
!     ------------------------------------------------------------------
    call reliem(' ', nomail, 'NO_MAILLE', motfac, iocc, &
                nbmc, motcle, typmcl, mesmai, nbma)
!
    if (nbma .eq. 0) then
        call utmess('F', 'ELEMENTS_66')
    end if
    call jeveuo(mesmai, 'L', jmail)
    call wkvect(mafour, 'V V I', nbma, jcour2)
    do im = 1, nbma
        ima = char8_to_int(zk8(jmail-1+im))
        zi(jcour2-1+im) = ima
    end do
!
!
!     ------------------------------------------------------------------
!     --- VERIFICATION DE L'EXISTENCE DES MAILLES
!     --- VERIFICATION QUE LES MAILLES SONT TOUTES SEG2, SEG3, SEG4 :
!     ------------------------------------------------------------------
    typmp = ' '
    ier = 0
    call jelira(typp, 'LONMAX', nbmato)
    do im = 1, nbma
        nomma = zk8(jmail-1+im)
        lexist = char8_to_int(nomma) .gt. nbmato
        if (lexist) then
            ier = ier+1
            call utmess('E', 'ELEMENTS5_19', sk=nomma, si=iocc)
        else
            ibid = char8_to_int(nomma)
            jtypm = iatyma-1+ibid
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypm)), typm)
            if (typm(1:3) .ne. 'SEG') then
                if (nomcmd .ne. 'DEFI_GROUP') then
                    call utmess('F', 'RUPTURE0_63')
                end if
                ier = ier+1
                call utmess('E', 'ELEMENTS5_20', sk=nomma, si=iocc)
            end if
            if (im .gt. 1) then
                if (typm .ne. typmp) then
                    ier = ier+1
                    call utmess('E', 'ELEMENTS5_21', sk=nomma, si=iocc)
                end if
            end if
            typmp = typm
        end if
    end do
    if (ier .gt. 0) then
        call utmess('F', 'ELEMENTS5_15', si=iocc)
    end if
!
!
! --- LECTURE DU NOM DU NOEUD ORIGINE (S'IL EST FOURNI)
    call getvtx(motfac, 'NOEUD_ORIG', iocc=iocc, nbval=0, nbret=n1)
    call getvtx(motfac, 'GROUP_NO_ORIG', iocc=iocc, nbval=0, nbret=n2)
    if (n1 .ne. 0) then
        call getvtx(motfac, 'NOEUD_ORIG', iocc=iocc, scal=ndorig, nbret=n1)
    else if (n2 .ne. 0) then
        call getvtx(motfac, 'GROUP_NO_ORIG', iocc=iocc, scal=nogrp, nbret=n2)
        call utnono(' ', nomail, 'NOEUD', nogrp, ndorig, &
                    iret)
        if (iret .eq. 10) then
            call utmess('F', 'ELEMENTS_67', sk=nogrp)
        else if (iret .eq. 1) then
            valk(1) = 'GROUP_NO_ORIG'
            valk(2) = ndorig
            call utmess('A', 'ELEMENTS5_17', nk=2, valk=valk)
        end if
    else
        ndorig = ' '
    end if
!
!
! --- LECTURE DU NOM DU NOEUD EXTREMITE (S'IL EST FOURNI)
    call getvtx(motfac, 'NOEUD_EXTR', iocc=iocc, nbval=0, nbret=n1)
    call getvtx(motfac, 'GROUP_NO_EXTR', iocc=iocc, nbval=0, nbret=n2)
    if (n1 .ne. 0) then
        call getvtx(motfac, 'NOEUD_EXTR', iocc=iocc, scal=ndextr, nbret=n1)
    else if (n2 .ne. 0) then
        call getvtx(motfac, 'GROUP_NO_EXTR', iocc=iocc, scal=nogrp, nbret=n2)
        call utnono(' ', nomail, 'NOEUD', nogrp, ndextr, &
                    iret)
        if (iret .eq. 10) then
            call utmess('F', 'ELEMENTS_67', sk=nogrp)
        else if (iret .eq. 1) then
            valk(1) = 'GROUP_NO_EXTR'
            valk(2) = ndextr
            call utmess('A', 'ELEMENTS5_17', nk=2, valk=valk)
        end if
    else
        ndextr = ' '
    end if
!
!
!     ------------------------------------------------------------------
!     --- CONSTRUCTION
!     --- 1 - VECTEUR DE TRAVAIL LOCAL CONTENANT LES NOEUDS EXTREMITES
!     ---     DE CHAQUE MAILLE
!     --- 2 - VECTEUR DE TRAVAIL LOCAL CONTENANT POUR CHAQUE NOEUD
!     ---     DU MAILLAGE :
!     ---     - 0 SI LE NOEUD N'APPARTIENT AUX MAILLES
!     ---     - 1 SI INTERNE (APPARTIENT A DEUX MAILLES)
!     ---     - 2 SI EXTREMITE
!     ------------------------------------------------------------------
    AS_ALLOCATE(vi=noeuds_extrem, size=2*nbma)
    AS_ALLOCATE(vi=type_noeud, size=nbnot)
    do im = 1, nbma
        call i2extf(zi(jcour2-1+im), 1, conec(1:15), typp(1:16), nig, &
                    nid)
        noeuds_extrem(im) = nig
        noeuds_extrem(nbma+im) = nid
        type_noeud(nig) = type_noeud(nig)+1
        type_noeud(nid) = type_noeud(nid)+1
    end do
!
!
! --- VERIFICATION QUE LA LIGNE EST CONTINUE ET UNIQUE
    n1 = 0
    n2 = 0
    bug = .false.
    do im = 1, nbnot
!        COMPTAGE DES EXTREMITES
        if (type_noeud(im) .eq. 1) n1 = n1+1
!        COMPTAGE NOEUDS APPARTENANT A PLUS DE DEUX MAILLES
        if (type_noeud(im) .gt. 2) n2 = n2+1
    end do
!     IL NE PEUT Y AVOIR QUE 2 NOEUDS EXTREMITES
    if (n1 .gt. 2) bug = .true.
!     IL NE DOIT PAS Y AVOIR DE NOEUDS APPARTENANT A PLUS DE DEUX
!     MAILLES
    if (n2 .ne. 0) bug = .true.
    if (bug) then
        call utmess('F', 'ELEMENTS5_16')
    end if
!
!
!     -- CALCUL DE VECORI:
!     --------------------
    vecori(1) = 0.d0
    vecori(2) = 0.d0
    vecori(3) = 0.d0
    if (nomcmd .eq. 'DEFI_GROUP' .and. motfac .eq. 'CREA_GROUP_NO') then
        call getvr8(motfac, 'VECT_ORIE', iocc=iocc, nbval=3, vect=vecori, &
                    nbret=n1)
        if ((ndorig .eq. ndextr) .and. (ndorig .ne. ' ')) then
            if (n1 .le. 0) then
                call utmess('A', 'ELEMENTS_70')
            end if
        end if
    end if
!
!
!
!   ------------------------------------------------------------------
!   --- verification du noeud extremite :
!   ------------------------------------------------------------------
    if (ndextr .ne. ' ') then
        nunori = char8_to_int(ndextr)
!
!       ON VERIFIE QU'IL S'AGIT BIEN D'UNE EXTREMITE
        trouv = 0
        do im = 1, nbma
            if (noeuds_extrem(im) .eq. nunori) trouv = trouv+1
            if (noeuds_extrem(nbma+im) .eq. nunori) trouv = trouv+1
        end do
!
        if (trouv .eq. 0) then
            call utmess('F', 'ELEMENTS_68', sk=ndorig)
        end if
        if (typlig .eq. 'FERME') then
            if (trouv .ne. 2) then
                call utmess('F', 'ELEMENTS_69', sk=ndextr)
            end if
        else
            if (.not. (typlig .eq. ' ' .and. ndorig .eq. ndextr)) then
                if (trouv .ne. 1) then
                    call utmess('F', 'ELEMENTS_69', sk=ndextr)
                end if
            end if
        end if
    end if
!
!
!
!   ------------------------------------------------------------------
!   --- Verification du noeud origine :
!   ------------------------------------------------------------------
    if (ndorig .ne. ' ') then
        nunori = char8_to_int(ndorig)
!
!       ON VERIFIE QU'IL S'AGIT BIEN D'UNE EXTREMITE
        trouv = 0
        do im = 1, nbma
            if (noeuds_extrem(im) .eq. nunori) trouv = trouv+1
            if (noeuds_extrem(nbma+im) .eq. nunori) trouv = trouv+1
        end do
!
        if (trouv .eq. 0) then
            call utmess('F', 'ELEMENTS_68', sk=ndorig)
        end if
        if (typlig .eq. 'FERME') then
            if (trouv .ne. 2) then
                call utmess('F', 'ELEMENTS_69', sk=ndorig)
            end if
        else
            if (.not. (typlig .eq. ' ' .and. ndorig .eq. ndextr)) then
                if (trouv .ne. 1) then
                    call utmess('F', 'ELEMENTS_69', sk=ndorig)
                end if
            end if
        end if
!
    else
!
!       ------------------------------------------------------------------
!       --- Si l'origine n'est pas donnee, on en cherche une :
!       ------------------------------------------------------------------
        AS_ALLOCATE(vi=noeud_apparies, size=2*nbma)
        noeud_apparies(:) = 0

!       -- parcours de l'ensemble des noeuds
        do in = 1, nbma*2
            if (noeud_apparies(in) .ne. 0) goto 80
            nunori = noeuds_extrem(in)

            do nd = in+1, nbma*2
                if (noeuds_extrem(nd) .eq. nunori) then
                    noeud_apparies(nd) = 1
                    goto 80
                end if
            end do

!           -- nunori n'apparait qu'une fois : c'est l'origine
            goto 100
!
80          continue
        end do

!       La ligne est peut etre fermee. On peut s'en sortir si ORIGINE='SANS' :
        erreur = .true.
        if (getexm(motfac, 'ORIGINE') .eq. 1) then
            call getvtx(motfac, 'ORIGINE', iocc=iocc, scal=orig, nbret=n3)
            if (n3 .eq. 1) then
                ASSERT(orig .eq. 'SANS')
                erreur = .false.
                nunori = noeuds_extrem(1)
!               -- On va verifier (grossierement) que la ligne est fermee.
!                  (on verifie seulement que chaque noeud est utilise 2 fois) :
                AS_ALLOCATE(vi=compteur, size=nbnot)
                compteur(:) = 0
                do im = 1, nbma
                    compteur(noeuds_extrem(im)) = compteur(noeuds_extrem(im))+1
                    compteur(noeuds_extrem(nbma+im)) = compteur(noeuds_extrem(nbma+im))+1
                end do
                do in = 1, nbnot
                    if (compteur(in) .eq. 0 .or. compteur(in) .eq. 2) cycle
                    erreur = .true.
                end do
                AS_DEALLOCATE(vi=compteur)
            end if
        end if
        if (erreur) call utmess('F', 'ELEMENTS_71')

100     continue
!
        ndorig = int_to_char8(nunori)
        call utmess('I', 'ELEMENTS_72', sk=ndorig)
!
        AS_DEALLOCATE(vi=noeud_apparies)
    end if

    AS_DEALLOCATE(vi=noeuds_extrem)
    AS_DEALLOCATE(vi=type_noeud)
    call jedetr('&&CGNOOR.MES_MAILLES')

    call jedema()
end subroutine
