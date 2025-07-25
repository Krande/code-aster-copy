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

subroutine irchme(ifichi, chanom, partie, nochmd, noresu, &
                  nomsym, typech, numord, nbrcmp, nomcmp, &
                  nbnoec, linoec, nbmaec, limaec, lvarie, &
                  sdcarm, carael, paraListNb, paraListName, &
                  nbCmpDyna, lfichUniq, codret)
!_______________________________________________________________________
!        IMPRESSION DU CHAMP CHANOM NOEUD/ELEMENT ENTIER/REEL
!        AU FORMAT MED
!     ENTREES:
!        IFICHI : UNITE LOGIQUE D'IMPRESSION DU CHAMP
!        CHANOM : NOM ASTER DU CHAM A ECRIRE
!        PARTIE : IMPRESSION DE LA PARTIE IMAGINAIRE OU REELLE POUR
!                  UN CHAMP COMPLEXE AU FORMAT CASTEM OU GMSH OU MED
!        NORESU : NOM DU RESULTAT D'OU PROVIENT LE CHAMP A IMPRIMER.
!        NOMSYM : NOM SYMBOLIQUE DU CHAMP
!        TYPECH : TYPE DU CHAMP
!        NUMORD : NUMERO D'ORDRE DU CHAMP DANS LE RESULTAT_COMPOSE.
!        NBRCMP : NOMBRE DE COMPOSANTES A ECRIRE
!        NOMCMP : NOMS DES COMPOSANTES A ECRIRE
!        NBNOEC : NOMBRE DE NOEUDS A ECRIRE (O, SI TOUS LES NOEUDS)
!        LINOEC : LISTE DES NOEUDS A ECRIRE SI EXTRAIT
!        NBMAEC : NOMBRE DE MAILLES A ECRIRE (0, SI TOUTES LES MAILLES)
!        LIMAEC : LISTE DES MAILLES A ECRIRE SI EXTRAIT
!        SDCARM : CARA_ELEM (UTILE POUR LES SOUS-POINTS)
! In  paraListNb       : length of list of parameter names
! Ptr paraListName     : pointer to the list of parameter names
!     SORTIES:
!        CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_______________________________________________________________________
!
!     ARBORESCENCE DE L'ECRITURE DES CHAMPS AU FORMAT MED :
!  IRCH19
!  IRCHME
!  MDNOCH RSADPA  IRCNME IRCEME
!                   .    .
!                    .  .
!                   IRCAME
!                    .  .
!                   .    .
!  MDNOMA MDEXMA IRMAIL UTLICM LRMTYP IRCMPR MDEXCH EFOUVR ...
!                   ... IRCMCC IRCMPG IRCMVA IRCMEC EFFERM
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/utflsh.h"
#include "asterfort/exisd.h"
#include "asterfort/copisd.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/irceme.h"
#include "asterfort/ircnme.h"
#include "asterfort/irmpav.h"
#include "asterfort/irmeta.h"
#include "asterfort/irvari.h"
#include "asterfort/detrsd.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnum.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: noresu, typech, sdcarm, carael
    character(len=16) :: nomsym
    character(len=19) :: chanom, ligrel
    character(len=24) :: nocelk
    character(len=*) :: nomcmp(*), partie
    integer(kind=8), intent(in) :: paraListNb
    character(len=8), pointer :: lgrf(:) => null()
    character(len=16), pointer :: paraListName(:)
    integer(kind=8) :: numord, nbrcmp, ifichi, iret
    integer(kind=8) :: nbnoec, nbmaec, icelk
    integer(kind=8) :: linoec(*), limaec(*)
!
    aster_logical :: lvarie, lfichUniq
!
    integer(kind=8), intent(inout) :: nbCmpDyna
    integer(kind=8) :: codret
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    integer(kind=8) :: ednono
    parameter(ednono=-1)
    integer(kind=8) :: ednopt
    parameter(ednopt=-1)
!
    integer(kind=8) :: ifm, nivinf, numpt, iaux, nbgrel, jmaille, j1, n1
    integer(kind=8) :: nbma, igr, iel, ite, ima, codret_vari
!
    character(len=8) :: saux08, modele, fauxmodele, ma
    character(len=64) :: nochmd
    real(kind=8) :: start_time, end_time
!
    real(kind=8) :: instan
!
!====
! 1. PREPARATIFS
!====
!
    call infniv(ifm, nivinf)
    codret = 0
    codret_vari = 0
    fauxmodele = ' '
!
100 format(/, 81('='), /, 81('='),/)
101 format(81('-'),/)
    if (nivinf .gt. 1) then
        call cpu_time(start_time)
        call utflsh(codret)
        write (ifm, 100)
        call utmess('I', 'MED_90', sk=chanom)
    end if
!
! 1.1. ==> NOM DU CHAMP DANS LE FICHIER MED
!
    saux08 = noresu
!
    if (codret .eq. 0) then
        if (nivinf .gt. 1) then
            write (ifm, 110) saux08
            write (ifm, 111) nomsym
        end if
        if (nivinf .gt. 1) then
            write (ifm, 113) typech
            write (ifm, 114) nochmd
        end if
    else
        call utmess('A', 'MED_91')
        call utmess('A', 'MED_44', sk=chanom)
        call utmess('A', 'MED_45', sk=noresu)
    end if
!
110 format(1x, 'RESULTAT           : ', a8)
111 format(1x, 'CHAMP              : ', a16)
113 format(1x, 'TYPE DE CHAMP      : ', a)
114 format(3x, '==> NOM MED DU CHAMP : ', a64,/)
!
! 1.2. ==> INSTANT CORRESPONDANT AU NUMERO D'ORDRE
!
    if (codret .eq. 0) then
!
        if (noresu .ne. ' ') then
            instan = 999.999d0
!         -- DANS UN EVOL_NOLI, IL PEUT EXISTER INST ET FREQ.
!            ON PREFERE INST :
            call jenonu(jexnom(noresu//'           .NOVA', 'INST'), iret)
            if (iret .ne. 0) then
                call rsadpa(noresu, 'L', 1, 'INST', numord, &
                            0, sjv=iaux, styp=saux08, istop=0)
                instan = zr(iaux)
            else
                call jenonu(jexnom(noresu//'           .NOVA', 'FREQ'), iret)
                if (iret .ne. 0) then
                    call rsadpa(noresu, 'L', 1, 'FREQ', numord, &
                                0, sjv=iaux, styp=saux08, istop=0)
                    instan = zr(iaux)
                else
                    call jenonu(jexnom(noresu//'           .NOVA', 'CHAR_CRIT'), iret)
                    if (iret .ne. 0) then
                        call rsadpa(noresu, 'L', 1, 'CHAR_CRIT', numord, &
                                    0, sjv=iaux, styp=saux08, istop=0)
                        instan = zr(iaux)
                    end if
                end if
            end if
            numpt = numord
            if (paraListNb .gt. 0) then
                call irmpav(noresu, ifichi, paraListNb, paraListName, numpt, numord, instan)
            end if
!
        else
!
            numord = ednono
            numpt = ednopt
!
        end if
!
    end if
!
! 1.3. ==> recherche du nom du modele (pour les champs ELGA) :
!
    if (codret .eq. 0) then
!
        if (typech(1:4) .ne. 'ELGA' .and. typech(1:4) .ne. 'ELNO') then
            modele = ' '
        else
            nocelk = chanom//'.CELK'
            call jeveuo(nocelk, 'L', icelk)
            ligrel = zk24(icelk) (1:19)
            call jeveuo(ligrel//'.LGRF', 'L', vk8=lgrf)
            modele = lgrf(2)

            call exisd('MODELE', modele, iret)
            if (iret .eq. 0) then
                if (noresu .ne. ' ') then
                    call rsadpa(noresu, 'L', 1, 'MODELE', numord, &
                                0, sjv=iaux, styp=saux08, istop=0)
                    modele = zk8(iaux)
                    call exisd('MODELE', modele, iret)
                end if
            end if
            if (iret .eq. 0) then
!            -- En absence d'un modele on va en construire un faux pour l'impression.
                fauxmodele = '&&IRCHME'
                call dismoi('NOM_MAILLA', chanom, 'CHAMP', repk=ma)
                call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
                call wkvect(fauxmodele//'.MAILLE', 'V V I', nbma, jmaille)
                call jelira(ligrel//'.LIEL', 'NUTIOC', nbgrel)
                do igr = 1, nbgrel
                    call jelira(jexnum(ligrel//'.LIEL', igr), 'LONMAX', n1)
                    call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', j1)
                    ite = zi(j1-1+n1)
                    do iel = 1, n1-1
                        ima = zi(j1-1+iel)
                        ASSERT(ima .ge. 0 .and. ima .le. nbma)
                        if (ima .gt. 0) zi(jmaille-1+ima) = ite
                    end do
                end do
                call copisd('LIGREL', 'V', ligrel, fauxmodele//'.MODELE')
                modele = fauxmodele
            end if
        end if
!
    end if
!
!====
! 2. ECRITURE DANS LE FICHIER MED
!====
!
    codret_vari = 0
    if (codret .eq. 0) then
!
        if (typech(1:4) .eq. 'NOEU') then
            call ircnme(ifichi, nochmd, chanom, typech, modele, &
                        nbrcmp, nomcmp, partie, numpt, instan, &
                        numord, nbnoec, linoec, sdcarm, carael, &
                        nomsym, lfichUniq, codret)
        else if (typech(1:2) .eq. 'EL') then
            if ((nomsym .eq. 'VARI_ELGA') .and. lvarie) then
                call irvari(ifichi, nochmd, chanom, typech, modele, &
                            nbrcmp, nomcmp, partie, numpt, instan, &
                            numord, nbmaec, limaec, noresu, sdcarm, &
                            carael, lfichUniq, codret_vari)
            end if
            if ((nomsym .eq. 'META_ELNO') .and. lvarie) then
                call irmeta(ifichi, nochmd, chanom, typech, modele, &
                            nbrcmp, nomcmp, partie, numpt, instan, &
                            numord, nbmaec, limaec, noresu, sdcarm, &
                            lfichUniq, codret_vari)
            end if
            call irceme(ifichi, nochmd, chanom, typech, modele, &
                        nbrcmp, nomcmp, ' ', partie, numpt, &
                        instan, numord, nbmaec, limaec, sdcarm, &
                        carael, nomsym, nbCmpDyna, lfichUniq, codret)
            if (codret_vari .ne. 0 .and. codret .eq. 0) then
                codret = codret_vari
            end if
        else if (typech(1:4) .eq. 'CART') then
            call irceme(ifichi, nochmd, chanom, typech, modele, &
                        nbrcmp, nomcmp, ' ', partie, numpt, &
                        instan, numord, nbmaec, limaec, sdcarm, &
                        carael, nomsym, nbCmpDyna, lfichUniq, codret)
        else
            codret = 1
            call utmess('A', 'MED_92', sk=typech(1:4))
        end if
!
    end if
!
!====
! 3. BILAN
!====
!
    if (codret .ne. 0 .and. &
        codret .ne. 100 .and. codret .ne. 200 .and. codret .ne. 300 .and. codret .ne. 400) then
        call utmess('A', 'MED_89', sk=nomsym)
    end if
!
    if (nivinf .gt. 1) then
        call cpu_time(end_time)
        call utmess('I', 'MED_93', sk=chanom, sr=end_time-start_time)
        write (ifm, 100)
        call utflsh(codret)
        write (ifm, 101)
    end if
!
    call detrsd('MODELE', fauxmodele)
!
end subroutine
