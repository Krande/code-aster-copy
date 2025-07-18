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
subroutine mdexcc(nofimd, idfimd, nochmd, nbcmpc, nomcmc, &
                  existc, nbcmfi, nmcmfi, codret)
! person_in_charge: nicolas.sellenet at edf.fr
!_____________________________________________________________________
!        FORMAT MED : EXISTENCE D'UN CHAMP - EST-IL CREE DANS UN FICHIER
!               - -   --             -              -
!_______________________________________________________________________
! .        .     .        .                                            .
! .  NOM   . E/S . TAILLE .           DESCRIPTION                      .
! .____________________________________________________________________.
! . NOFIMD .  E  .   1    . NOM DU FICHIER MED                         .
! . NOFIMD .  E  .   1    . OU NUMERO DU FICHIER DEJA OUVERT           .
! . NOCHMD .  E  .   1    . NOM DU CHAMP MED VOULU                     .
! . NBCMPC .  E  .   1    . NOMBRE DE COMPOSANTES A CONTROLER          .
!          .     .        . S'IL EST NUL, ON NE CONTROLE RIEN          .
! . NOMCMC .  E  .   *    . SD DES NOMS DES COMPOSANTES A CONTROLER    .
! . NOMAMD .  E  .   1    . NOM DU MAILLAGE MED ASSOCIE                .
! . EXISTC .  S  .   1    . 0 : LE CHAMP N'EST PAS CREE                .
! .        .     .        . >0 : LE CHAMP EST CREE AVEC :              .
! .        .     .        . 1 : LES COMPOSANTES VOULUES SONT PRESENTES .
! .        .     .        . 2 : LES COMPOSANTES VOULUES NE SONT PAS    .
! .        .     .        .     TOUTES ENREGISTREES                    .
! . NBCMFI .  S  .   1    . NOMBRE DE COMPOSANTES DANS LE FICHIER      .
! . NMCMFI .  S  .   1    . SD DU NOM DES COMPOSANTES DANS LE FICHIER  .
! . CODRET .  S  .    1   . CODE DE RETOUR DES MODULES                 .
! ______________________________________________________________________
!
!====
! 0. DECLARATIONS ET DIMENSIONNEMENT
!====
!
    use as_med_module, only: as_med_open
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_mfdfdi.h"
#include "asterfort/as_mfdnfc.h"
#include "asterfort/as_mfdnfd.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/codent.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: nofimd, nochmd
    character(len=*) :: nomcmc, nmcmfi
!
    integer(kind=8) :: nbcmpc, existc, nbcmfi
!
    integer(kind=8) :: codret
!
! 0.2. ==> COMMUNS
!
!
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='MDEXCC')
!
    integer(kind=8) :: edlect
    integer(kind=8) :: vali(2)
    parameter(edlect=0)
    integer(kind=8) :: mfloat
    parameter(mfloat=6)
!
!
    integer(kind=8) :: lnochm, jnbcha, jnocha, jcmpch
    med_idt :: idfimd
    integer(kind=8) :: nbcham, nbcha2
    integer(kind=8) :: iaux, jaux, kaux, iret, iouv
    integer(kind=8) :: adncmp, aducmp, adncmc, adncfi, nseqca
    aster_logical :: ficexi, dejouv
!
    character(len=8) :: saux08
    character(len=16) :: saux16
    character(len=24) :: nonbch, nonoch, nocmch, nomcmp
    character(len=64) :: saux64
! ______________________________________________________________________
!
    existc = 0
    nbcmfi = -1
    codret = 0
!
    nonbch = '&&'//nompro//'.NB_CHAMPS'
    nonoch = '&&'//nompro//'.NOM_CHAMPS'
    nocmch = '&&'//nompro//'.CMP_CHAMPS'
!
!====
! 1. ON CREE LES TABLEAUX CONTENANT LE NOMBRE DE CHAMPS, LE NOM
!     DES CHAMPS ET SES COMPOSANTES
!====
!
! 1.1 ==> OUVERTURE DU FICHIER S'IL N'EST PAS DEJA OUVERT
!
    inquire (file=nofimd, exist=ficexi)
!
    if (.not. ficexi) goto 999
!
    if (idfimd .eq. 0) then
        call as_med_open(idfimd, nofimd, edlect, iouv)
        dejouv = .false.
    else
        dejouv = .true.
        iouv = 0
    end if
!
! 1.2 ==> SI ON A PU OUVRIR LE FICHIER, ON COMMENCER A LE LIRE
!
    if (iouv .eq. 0) then
!
! 1.3 ==> LECTURE DU NOMBRE DE CHAMPS
!
        call as_mfdnfd(idfimd, nbcham, codret)
        if (codret .ne. 0) then
            saux08 = 'mfdnfd'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
        if (nbcham .le. 0) goto 30
!
        call jeexin(nonbch, iret)
        if (iret .ne. 0) then
!
!         SI LE VECTEUR EXISTE MAIS QUE LE NOMBRE DE
!         CHAMP A CHANGE, IL FAUT RE CONSTRUIRE CES
!         VECTEURS
            call jeveuo(nonbch, 'L', jnbcha)
            nbcha2 = zi(jnbcha)
            if (nbcha2 .ne. nbcham) then
                call jeveuo(nocmch, 'L', jcmpch)
                do iaux = 1, nbcha2
                    nomcmp = zk24(jcmpch+iaux-1)
                    call jedetr(nomcmp)
                end do
                call jedetr(nonbch)
                call jedetr(nonoch)
                call jedetr(nocmch)
                iret = 0
            end if
        end if
!
        if (iret .eq. 0) then
!
! 1.4 ==> ALLOCATION DU TABLEAU CONTENANT : - LE NOMBRE DE CHAMPS
!                                           - LE NOM DES CHAMPS
!                                           - LES COMPOSANTES
!
            call wkvect(nonbch, 'V V I', 1, jnbcha)
            call wkvect(nonoch, 'V V K80', nbcham, jnocha)
            call wkvect(nocmch, 'V V K24', nbcham, jcmpch)
            zi(jnbcha) = nbcham
!
! 1.5 ==> POUR CHAQUE CHAMP ON CHERCHE SON NOM ET SES COMPOSANTES
!
            do 10, iaux = 1, nbcham
                call as_mfdnfc(idfimd, iaux, nbcmfi, codret)
                if (codret .ne. 0) then
                    saux08 = 'mfdnfc'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                call codent(iaux, 'G', saux08)
                call wkvect('&&'//nompro//saux08//'N', 'V V K16', nbcmfi, adncmp)
                call wkvect('&&'//nompro//saux08//'U', 'V V K16', nbcmfi, aducmp)
                saux64 = ' '
!
! 1.5.1 ==> LECTURE DU NOM DU CHAMP MED ET DE SES COMPOSANTES
!
                call as_mfdfdi(idfimd, iaux, saux64, jaux, zk16(adncmp), &
                               zk16(aducmp), nseqca, codret)
                if (codret .ne. 0 .or. jaux .ne. mfloat) then
                    vali(1) = iaux
                    if (codret .ne. 0) then
                        saux08 = 'mfdfdi'
                        call utmess('F', 'DVP_97', sk=saux08, si=codret)
                    end if
                    if (saux64 .eq. nochmd .and. jaux .ne. mfloat) then
                        vali(1) = jaux
                        call utmess('A+', 'MED_84', si=vali(1))
                        call utmess('F', 'MED_75')
                    end if
                end if
!
! 1.5.2 ==> RECOPIE DANS LES CHAMPS ALLOUES
!
                zk80(jnocha+iaux-1) = saux64
                zk24(jcmpch+iaux-1) = '&&'//nompro//saux08//'N'
                call jedetr('&&'//nompro//saux08//'U')
10              continue
                end if
!
!====
! 2. LE CHAMP EST-IL PRESENT ?
!====
!
! 2.1. ==> NBCHAM : NOMBRE DE CHAMPS DANS LE FICHIER
!
                call jeveuo(nonoch, 'L', jnocha)
                call jeveuo(nocmch, 'L', jcmpch)
!
! 2.2. ==> RECHERCHE DU CHAMP VOULU
!
                lnochm = lxlgut(nochmd)
!
                do 22, iaux = 1, nbcham
!
                    saux64 = zk80(jnocha+iaux-1)
                    nomcmp = zk24(jcmpch+iaux-1)
!
! 2.2.3. ==> COMPARAISON DU NOM DU CHAMP
!
                    jaux = lxlgut(saux64)
!
                    if (jaux .eq. lnochm) then
                        if (saux64(1:jaux) .eq. nochmd(1:lnochm)) then
                            existc = 1
                        end if
                    end if
!
! 2.2.4. ==> C'EST LE BON CHAMP. CONTROLE DU NOM DES COMPOSANTES
!
                    if (existc .eq. 1) then
!
! 2.2.4.1. ==> TRANSFERT DES NOMS DES COMPOSANTES DANS LE TABLEAU
!              DE SORTIE
!
                        call jelira(nomcmp, 'LONMAX', nbcmfi)
                        call jeveuo(nomcmp, 'L', adncmp)
!
                        call wkvect(nmcmfi, 'V V K16', nbcmfi, adncfi)
!
                        do 2241, kaux = 0, nbcmfi-1
                            zk16(adncfi+kaux) = zk16(adncmp+kaux)
2241                        continue
!
! 2.2.4.2. ==> TEST DES NOMS DES COMPOSANTES
!
                            if (nbcmpc .gt. nbcmfi) then
!
                                existc = 2
!
                            else if (nbcmpc .gt. 0) then
!
                                call jeveuo(nomcmc, 'L', adncmc)
!
!      ZK16(ADNCMC+JAUX) : NOM DE LA (JAUX+1)-EME COMPOSANTE A CONTROLER
!      ZK16(ADNCMP+KAUX) : NOM DE LA (KAUX+1)-EME COMPOSANTE DU CHAMP
!
                                do 2242, jaux = 0, nbcmpc-1
!
                                    saux16 = zk16(adncmc+jaux)
!
                                    do 2243, kaux = 0, nbcmfi-1
                                    if (saux16 .eq. zk16(adncmp+kaux)) then
                                        goto 2242
                                    end if
2243                                continue
!
!           AUCUNE COMPOSANTE DU CHAMP LU NE CORRESPOND A LA COMPOSANTE
!           SOUHAITEE
!
                                    existc = 2
!
2242                                continue
!
                                    end if
!
                                    goto 30
!
                                    end if
!
22                                  continue
!
!====
! 3. LA FIN
!====
!
30                                  continue
!
!====
! 4. ==> FERMETURE DU FICHIER S'IL Y A BESOIN
!====
!
                                    if (.not. dejouv) then
                                        call as_mficlo(idfimd, codret)
                                        if (codret .ne. 0) then
                                            saux08 = 'mficlo'
                                            call utmess('F', 'DVP_97', sk=saux08, si=codret)
                                        end if
                                        idfimd = 0
                                    end if
                                    end if
!
999                                 continue
!
                                    end subroutine
