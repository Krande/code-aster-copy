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
subroutine lrcmva(ntvale, nbvato, ntproa, lgproa, ncmprf, &
                  nomcmr, nbcmfi, nmcmfi, nbcmpv, ncmpvm, &
                  numcmp, nochmd, adsl, adsv, codret)
!_____________________________________________________________________
! person_in_charge: nicolas.sellenet at edf.fr
! ======================================================================
!     LECTURE D'UN CHAMP - FORMAT MED - CREATION DES VALEURS AUX NOEUDS
!     -    -       -              -                  --
!-----------------------------------------------------------------------
!     ENTREES:
!       NTVALE : TABLEAU QUI CONTIENT LES VALEURS LUES
!       NBVATO : NOMBRE DE VALEURS TOTAL
!       NTPROA : TABLEAU QUI CONTIENT LE PROFIL ASTER
!       LGPROA : LONGUEUR DU PROFIL ASTER
!                SI NUL, PAS DE PROFIL
!       NCMPRF : NOMBRE DE COMPOSANTES DU CHAMP DE REFERENCE
!       NOMCMR : NOMS DES COMPOSANTES DE REFERENCE
!       NBCMFI : NOMBRE DE COMPOSANTES DANS LE FICHIER      .
!       NMCMFI : NOM DES COMPOSANTES DANS LE FICHIER
!       NBCMPV : NOMBRE DE COMPOSANTES VOULUES
!                SI NUL, ON LIT LES COMPOSANTES A NOM IDENTIQUE
!       NCMPVM : LISTE DES COMPOSANTES VOULUES DANS MED
!       NUMCMP : TABLEAU DES NUMEROS DES COMPOSANTES VALIDES
!       NOCHMD : NOM MED DU CHAMP A LIRE
!     SORTIES:
!       ADSL, ADSV : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
!        CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!
!   REMARQUE :
!    LE TABLEAU DE VALEURS EST UTILISE AINSI : TV(NBCMFI,NVALUT)
!    EN FORTRAN, CELA CORRESPOND AU STOCKAGE MEMOIRE SUIVANT :
!    TV(1,1), TV(2,1), ..., TV(NBCMFI,1), TV(1,2), TV(2,2), ...,
!    TV(NBCMFI,2) , TV(1,3), TV(2,3), ..., TV(1,NVALUT), TV(2,NVALUT),
!    TV(NBCMFI,NVALUT)
!    C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!_____________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbvato, lgproa
    integer(kind=8) :: ncmprf, nbcmpv
    integer(kind=8) :: adsl, adsv
    integer(kind=8) :: codret
!
    character(len=*) :: nochmd
    character(len=*) :: nomcmr(*)
    character(len=*) :: ntvale, ntproa, ncmpvm, nmcmfi, numcmp
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='LRCMVA')
!
!
    integer(kind=8) :: iaux, jaux, kaux, laux, maux
    integer(kind=8) :: nrcmp, ncmpdb
    integer(kind=8) :: nuval
    integer(kind=8) :: nbcmfi
    integer(kind=8) :: adremp, advale, adncfi, adnucm, adncvm, adproa
    integer(kind=8) :: ifm, nivinf
!
    character(len=8) :: saux08
    character(len=24) :: ntcmpl
    character(len=24) :: valk(2)
!
    call infniv(ifm, nivinf)
!
!====
! 1. ON BOUCLE SUR LES NBCMFI COMPOSANTES DU CHAMP QUI A ETE LU
!    DANS LE FICHIER
!====
!
    codret = 0
!
!               12   345678   9012345678901234
    ntcmpl = '&&'//nompro//'.LOGIQUECMP     '
    call wkvect(ntcmpl, 'V V L', ncmprf, adremp)
!
    call jeveuo(ntvale, 'L', advale)
    if (nbcmpv .ne. 0) then
        call jeveuo(ncmpvm, 'L', adncvm)
    end if
    call jeveuo(nmcmfi, 'L', adncfi)
    call jeveuo(numcmp, 'L', adnucm)
!
    if (lgproa .ne. 0) then
        call jeveuo(ntproa, 'L', adproa)
    end if
!
    do iaux = 1, nbcmfi
!
        ncmpdb = 1
!
! 1.1. ==> REPERAGE DU NUMERO DE LA COMPOSANTE DU CHAMP ASTER DANS
!          LAQUELLE AURA LIEU LE TRANSFERT DE LA IAUX-EME COMPOSANTE LUE
!
110     continue
!
        nrcmp = 0
!
! 1.1.1. ==> QUAND ON VEUT UNE SERIE DE COMPOSANTES :
!            BOUCLE 111 : ON CHERCHE, DANS LA LISTE VOULUE EN MED,
!            NCMPVM, A QUEL ENDROIT SE TROUVE LA COMPOSANTE LUE EN
!            COURS DE TRAITEMENT, ZK16(ADNCFI-1+IAUX).
!            ON EXPLORE LA LISTE DE NCMPDB, QUI VAUT 1 AU DEBUT, A
!            NBCMPV, QUI EST LE NOMBRE DE COMPOSANTES VOULUES.
!
!            DEUX CAS DE FIGURE :
!            . SI AUCUN ELEMENT DE CETTE LISTE NE CORRESPOND A LA
!              COMPOSANTE LUE, CELA VEUT DIRE QUE CETTE COMPOSANTE N'EST
!              PAS SOUHAITEE. EN QUELQUE SORTE, ELLE A ETE LUE POUR
!              RIEN. ON PASSE A LA COMPOSANTE LUE SUIVANTE (cycle)
!
!            . QUAND ON TOMBE SUR UNE COMPOSANTE MED VOULUE IDENTIQUE,
!              ON VA REMPLIR LA COMPOSANTE ASTER ASSOCIEE AVEC LES
!              VALEURS DE LA COMPOSANTE LUE. POUR CELA :
!              . ON DEDUIT LE NUMERO DANS LA LISTE OFFICIELLE DE LA
!                COMPOSANTE ASTER ASSOCIEE, ZI(ADNUCM-1+JAUX).
!              . ON MEMORISE QUE LA COMPOSANTE ASTER EST REMPLIE GRACE
!                AU BOOLEN ZL(ADREMP+NRCMP-1).
!              . ON MEMORISE A QUEL ENDROIT DE LA LISTE VOULUE EN MED ON
!                EN EST. EN EFFET, LES VALEURS DE LA COMPOSANTE LUE
!                PEUVENT ETRE MISES DANS PLUSIEURS COMPOSANTES ASTER.
!                UNE FOIS LE TRANSFERT EFFECTUE, ON REPRENDRA
!                L'EXPLORATION DE LA LISTE A L'ENDROIT OU ON S'ETAIT
!                ARRETE : NCMPDB.
!              . ON FAIT LE TRANSFERT DES VALEURS (GOTO 12).
!
        if (nbcmpv .ne. 0) then
!
            do jaux = ncmpdb, nbcmpv
                if (zk16(adncvm-1+jaux) .eq. zk16(adncfi-1+iaux)) then
                    nrcmp = zi(adnucm-1+jaux)
                    ncmpdb = jaux+1
                    goto 12
                end if
            end do
            cycle
12          continue
!
! 1.1.2. ==> QUAND ON STOCKE A L'IDENTIQUE, ON RECHERCHE LE NUMERO DE
!            COMPOSANTE DE REFERENCE QUI A LE MEME NOM QUE LA
!            COMPOSANTE LUE. ON EST OBLIGE DE FAIRE CETTE RECHERCHE CAR
!            RIEN NE GARANTIT QUE L'ORDRE DES COMPOSANTES SOIT LE MEME
!            DANS LA REFERENCE ASTER ET DANS LE CHAMP ECRIT DANS LE
!            FICHIER.
        else
!
            saux08 = zk16(adncfi-1+iaux) (1:8)
            if (lxlgut(zk16(adncfi-1+iaux)) .gt. 8) then
                valk(1) = zk16(adncfi-1+iaux)
                valk(2) = saux08
                call utmess('A', 'MED_72', nk=2, valk=valk)
            end if
            nrcmp = indik8(nomcmr, saux08, 1, ncmprf)
!
        end if
!
! 1.2. ==> SI AUCUNE COMPOSANTE N'A ETE TROUVEE, MALAISE ...
!
        if (nrcmp .eq. 0) then
            call utmess('F', 'MED_73', sk=zk16(adncfi-1+iaux))
        end if
!
! 1.3. ==> TRANSFERT DES VALEURS DANS LA COMPOSANTE NRCMP
!
        zl(adremp+nrcmp-1) = .true.
!
! 1.3.1. ==> SANS PROFIL : ON PARCOURT LES NOEUDS DU PREMIER AU DERNIER
!
        if (lgproa .eq. 0) then
            kaux = -ncmprf+nrcmp-1
            laux = advale-nbcmfi+iaux-1
            do jaux = 1, nbvato
                kaux = kaux+ncmprf
                laux = laux+nbcmfi
                zl(adsl+kaux) = .true.
                zr(adsv+kaux) = zr(laux)
            end do
!
        else
!
! 1.3.2. ==> AVEC PROFIL : ON PARCOURT LES NOEUDS STOCKES DANS LE PROFIL
!            ON RECUPERE LEURS NUMEROS ET ON PLACE LA VALEUR AU BON
!            ENDROIT. ATTENTION, IL N'Y A AUCUNE RAISON POUR QUE LES
!            NUMEROS APPARAISENT DANS L'ORDRE, DONC IL FAUT RECALCULER
!            LA POSITION A CHAQUE FOIS.
!
            maux = -ncmprf+nrcmp-1
            laux = advale-nbcmfi+iaux-1
            do nuval = 0, lgproa-1
                jaux = zi(adproa+nuval)
                kaux = maux+jaux*ncmprf
                laux = laux+nbcmfi
                zl(adsl+kaux) = .true.
                zr(adsv+kaux) = zr(laux)
            end do
!
        end if
!
! 1.4. ==> QUAND ON VEUT UNE SERIE DE COMPOSANTES, ON REPREND
!          L'EXPLORATION DE LA LISTE VOULUE
!
        if (nbcmpv .ne. 0) then
            goto 110
        end if
!
    end do
!
!====
! 2. ON INFORME SUR LES COMPOSANTES QUI ONT ETE REMPLIES
!====
!
    kaux = 0
    do jaux = 1, ncmprf
        if (.not. zl(adremp+jaux-1)) then
            kaux = kaux+1
        end if
    end do
!
    if (kaux .gt. 0 .and. nivinf .gt. 1) then
        write (ifm, 2001) nochmd
        do jaux = 1, ncmprf
            if (zl(adremp+jaux-1)) then
                saux08 = nomcmr(jaux)
                write (ifm, 2002) saux08
            end if
        end do
        write (ifm, *) ' '
    end if
!
2001 format('CHAMP ', a)
2002 format('. LA COMPOSANTE LUE : ', a8, '.')
!
!====
! 3. MENAGE
!====
!
    call jedetr(ntcmpl)
!
end subroutine
