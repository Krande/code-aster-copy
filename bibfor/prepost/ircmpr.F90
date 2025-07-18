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
! aslint: disable=W1504
!
subroutine ircmpr(nofimd, typech, nbimpr, ncaimi, ncaimk, &
                  ncmprf, ncmpve, ntlcmp, nbvato, nbenec, &
                  lienec, adsd, adsl, nomaas, modele, &
                  typgeo, nomtyp, ntproa, chanom, sdcarm, &
                  field_type, nosdfu)
!_______________________________________________________________________
!     ECRITURE D'UN CHAMP -  FORMAT MED - CREATION DU PROFIL
!        -  -       -               -                 --
!_______________________________________________________________________
!     ENTREES :
!       NOFIMD : NOM DU FICHIER MED
!       TYPECH : TYPE DU CHAMP ('NOEU', 'ELNO', 'ELGA')
!       NCMPRF : NOMBRE DE COMPOSANTES DU CHAMP DE REFERENCE
!       NCMPVE : NOMBRE DE COMPOSANTES VALIDES EN ECRITURE
!       NTLCMP : SD DES NUMEROS DES COMPOSANTES VALIDES
!       NBVATO : NOMBRE DE VALEURS TOTALES
!       NBENEC : NOMBRE D'ENTITES A ECRIRE (O, SI TOUTES)
!       LIENEC : LISTE DES ENTITES A ECRIRE SI EXTRAIT
!       ADSK, D, ... : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
!       NOMAAS : SD MAILLAGE ASTER
!       MODELE : SD MODELE
!       TYPGEO : TYPE GEOMETRIQUE DE MAILLE ASSOCIEE AU TYPE ASTER
!       NOMTYP : NOM DES TYPES DE MAILLES ASTER
! In  field_type       : type of field (symbolic name in result datastructure)
!     SORTIES :
!       NBIMPR : NOMBRE D'IMPRESSIONS A REALISER
!       NCAIMI : STRUCTURE ASSOCIEE AU TABLEAU CAIMPI
!         CAIMPI : ENTIERS POUR CHAQUE IMPRESSION
!                  CAIMPI(1,I) = TYPE D'EF / MAILLE ASTER (0, SI NOEUD)
!                  CAIMPI(2,I) = NOMBRE DE POINTS (GAUSS OU NOEUDS)
!                  CAIMPI(3,I) = NOMBRE DE SOUS-POINTS
!                  CAIMPI(4,I) = NOMBRE DE COUCHES
!                  CAIMPI(5,I) = NOMBRE DE SECTEURS
!                  CAIMPI(6,I) = NOMBRE DE FIBTRES
!                  CAIMPI(7,I) = NOMBRE DE MAILLES A ECRIRE
!                  CAIMPI(8,I) = TYPE DE MAILLES ASTER (0, SI NOEUD)
!                  CAIMPI(9,I) = TYPE GEOMETRIQUE AU SENS MED
!                  CAIMPI(10,I) = NOMBRE TOTAL DE MAILLES IDENTIQUES
!       NCAIMK : STRUCTURE ASSOCIEE AU TABLEAU CAIMPK
!         CAIMPK : CARACTERES POUR CHAQUE IMPRESSION
!                  CAIMPK(1,I) = NOM DE LA LOCALISATION ASSOCIEE
!                  CAIMPK(2,I) = NOM DU PROFIL AU SENS MED
!       NTPROA : SD DU PROFIL ASTER. C'EST LA LISTE DES NUMEROS ASTER
!                DES NOEUDS OU DES ELEMENTS POUR LESQUELS LE CHAMP
!                EST DEFINI
!_______________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterc/utflsh.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/ircmpe.h"
#include "asterfort/ircmpn.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nbvato, ncmprf, ncmpve
    integer(kind=8) :: nbenec, adtyp2
    integer(kind=8) :: lienec(*)
    integer(kind=8) :: adsd, adsl
    integer(kind=8) :: nbimpr
    integer(kind=8) :: typgeo(*)
!
    character(len=*) :: nofimd
    character(len=*) :: ntlcmp, ntproa
    character(len=8) :: nomaas, modele, typech, sdcarm
    character(len=8) :: nomtyp(*), nosdfu
    character(len=19) :: chanom
    character(len=24) :: ncaimi, ncaimk
    character(len=16), intent(in) :: field_type
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='IRCMPR')
!
    integer(kind=8) :: ifm, nivinf, i, j, jco
    integer(kind=8) :: iaux, ima, nbno, nbma
    integer(kind=8) :: nbmail, iadcnx, ilcnx
    integer(kind=8) :: codret
    integer(kind=8) ::  adefma
    integer(kind=8) :: adcaii, adcaik
    integer(kind=8) :: adproa, adprom, adexic, adpror
    integer(kind=8) :: adnucm
    integer(kind=8) :: adauxi
    real(kind=8) :: start_time, end_time
!
    character(len=24) :: ntprom, exicmp, ntpror
    character(len=24) :: ntauxi
    character(len=80) :: caimpk(3)
    integer(kind=8), pointer :: noeu_centr(:) => null()
    integer(kind=8), pointer :: nadtypm(:) => null()
    integer(kind=8), pointer :: dtyp(:) => null()
!
!====
! 1. PREALABLES
!====
!
! 1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
!          -----------------------------------
!
    call infniv(ifm, nivinf)
!
    if (nivinf .gt. 1) then
        call cpu_time(start_time)
        write (ifm, 101) 'DEBUT DE '//nompro
        call utflsh(codret)
    end if
101 format(/, 4x, 10('='), a, 10('='),/)
!
!               12   345678   9012345678901234
    ntprom = '&&'//nompro//'.PROFIL_MED     '
    exicmp = '&&'//nompro//'.EXICMP         '
    ntauxi = '&&'//nompro//'.NUME_RECIPROQUE'
    ntpror = '&&'//nompro//'.PROFIL_RECIPROQ'
!
! 1.2. ==> ALLOCATIONS DES TABLEAUX DE RENUMEROTATIONS
!
    call wkvect(ntproa, 'V V I', nbvato, adproa)
    call wkvect(ntprom, 'V V I', nbvato, adprom)
    call wkvect(exicmp, 'V V L', nbvato, adexic)
!
    call jeveuo(ntlcmp, 'L', adnucm)
!
! 1.3. ==> COMPLEMENTS
!
! 1.3.1. ==> COMPLEMENTS POUR UN CHAMP AUX NOEUDS
!
    if (typech(1:4) .eq. 'NOEU') then
!
        nbimpr = 1
        iaux = 10*nbimpr
        call wkvect(ncaimi, 'V V I', iaux, adcaii)
        iaux = 3*nbimpr
        call wkvect(ncaimk, 'V V K80', iaux, adcaik)
!
!       ON CREE UN TABLEAU QUI PERMET DE DETECTER L'EXISTENCE DE NOEUDS
!       CENTRE (APPARTENANT AUX MAILLES DE TYPE TRIA7,QUAD9,PENTA18 OU
!       HEXA27)
!
        call jeveuo(nomaas//'.TYPMAIL', 'L', vi=dtyp)
        call jeveuo(nomaas//'.CONNEX', 'L', iadcnx)
        call jeveuo(jexatr(nomaas//'.CONNEX', 'LONCUM'), 'L', ilcnx)
        call dismoi('NB_NO_MAILLA', nomaas, 'MAILLAGE', repi=nbno)
        call dismoi('NB_MA_MAILLA', nomaas, 'MAILLAGE', repi=nbma)
        AS_ALLOCATE(vi=noeu_centr, size=nbno)
        do i = 1, nbno
            noeu_centr(i) = 0
        end do
!
        do i = 1, nbma
            if (dtyp(i) .eq. MT_TETRA15) then
                jco = iadcnx+zi(ilcnx+i-1)-1
                do j = 1, 5
                    noeu_centr(1+zi(jco+10+j-1)-1) = 1
                end do
            elseif (dtyp(i) .eq. MT_PYRAM19) then
                jco = iadcnx+zi(ilcnx+i-1)-1
                do j = 1, 6
                    noeu_centr(1+zi(jco+13+j-1)-1) = 1
                end do
            elseif (dtyp(i) .eq. MT_PENTA21) then
                jco = iadcnx+zi(ilcnx+i-1)-1
                do j = 1, 3
                    noeu_centr(1+zi(jco+18+j-1)-1) = 1
                end do
            end if
        end do
!
! 1.3.2. ==> COMPLEMENTS POUR DES CHAMPS AUX ELEMENTS
!
    else if (typech(1:2) .eq. 'EL') then
!
        call jeveuo(nomaas//'.TYPMAIL', 'L', vi=nadtypm)
        call jelira(nomaas//'.TYPMAIL', 'LONMAX', nbmail)
        call wkvect('&&IRCMPR.TYPMA', 'V V I', nbmail, adtyp2)
        do ima = 1, nbmail
            if (nadtypm(ima) .eq. MT_TETRA15) then
                zi(adtyp2+ima-1) = MT_TETRA10
            elseif (nadtypm(ima) .eq. MT_PYRAM19) then
                zi(adtyp2+ima-1) = MT_PYRAM13
            elseif (nadtypm(ima) .eq. MT_PENTA21) then
                zi(adtyp2+ima-1) = MT_PENTA18
            elseif (nadtypm(ima) .eq. MT_HEXA9) then
                zi(adtyp2+ima-1) = MT_HEXA8
            else
                zi(adtyp2+ima-1) = nadtypm(ima)
            end if
        end do
        if (typech(1:4) .eq. 'ELGA') then
            call jeveuo(modele//'.MAILLE', 'L', adefma)
        end if
        call wkvect(ntpror, 'V V I', nbvato, adpror)
        call wkvect(ntauxi, 'V V I', nbvato, adauxi)
!
! 1.3.3. ==> ERREUR
!
    else
!
        call utmess('F', 'MED_46', sk=typech)
!
    end if
!
!====
! 2. APPELS DES PROGRAMMES SPECIFIQUES
!====
!
    if (typech(1:4) .eq. 'NOEU') then
!
! 2.1. ==> LES NOEUDS
!
        call ircmpn(nofimd, ncmprf, ncmpve, zi(adnucm), zl(adexic), &
                    nbvato, nbenec, lienec, adsl, zi(adcaii), &
                    caimpk, zi(adproa), noeu_centr, nosdfu)
        zk80(adcaik:adcaik+2) = caimpk(1:3)
!
    else
!
! 2.2. ==> LES ELEMENTS
!
        if (typech(1:4) .eq. 'ELGA') then
            iaux = adefma
        else
            iaux = adtyp2
        end if
!
        call ircmpe(nofimd, ncmpve, zi(adnucm), zl(adexic), nbvato, &
                    nbenec, lienec, adsd, adsl, nbimpr, &
                    ncaimi, ncaimk, zi(iaux), zi(adtyp2), typgeo, &
                    nomtyp, typech, zi(adproa), zi(adprom), zi(adpror), &
                    zi(adauxi), chanom, sdcarm, field_type, nosdfu)
!
    end if
!
!====
! 3. LA FIN
!====
!
! --- MENAGE
    call jedetr(ntpror)
    call jedetr(ntprom)
    call jedetr(exicmp)
    call jedetr(ntauxi)
    AS_DEALLOCATE(vi=noeu_centr)
    call jedetr('&&IRCMPR.TYPMA')
!
    if (nivinf .gt. 1) then
        call cpu_time(end_time)
        write (ifm, *) '    ==========FIN DE '//nompro, " EN ", &
            end_time-start_time, " sec=========="
    end if
!
end subroutine
