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

subroutine lrvcpg(idfimd, nbpgm, nbpga, nomtm, typgeo, &
                  elrefa, fapg, nloc, locnam, permu, &
                  nutyma, nbsp, codret)
!
! person_in_charge: nicolas.sellenet at edf.fr
!     LECTURE FICHIER MED - VERIFICATION ET COMPARAISON DES PG ASTER/MED
!     -    -                -               -               --
!-----------------------------------------------------------------------
!
!     ROUTINE APPELEE PAR: LRMPGA
!
!     IN :
!       IDFIMD : IDENTIFIANT DU FICHIER MED
!       NBPGM  : NOMBRE DE PG MED
!       NBPGA  : NOMBRE DE PG ASTER
!       ELREFA : NOM DE L'ELEMENT DE REFERENCE ASTER
!       ELREFM : NOM DE L'ELEMENT DE REFERENCE MED
!       FAPG   : FAMILLE DE PG GLOBALE
!       NLOC   : NOMBRE DE LOCALISATIONS PRESENTES DANS LE FICHIER MED
!    IN/OUT:
!       PERMU  : TABLEAU (EVENTUEL) DES PERMUTATIONS DES PG
!                PERMU(I_PG_MED)=I_PG_ASTER
!    OUT :
!       NUTYMA : NUMERO DU TYPE DE MAILLE DE ELREFA
!       CODRET : CORRESPONDANCE DES PG ASTER/MED
!                 CODRET=0 -->  OK
!                 CODRET=1 -->  NECESSITE DES PERMUTATIONS
!                 CODRET=2 -->  NOOK (SOIT ABSENCE DE LOCALISATION, SOIT
!                                     AUCUNE CORRESPONDANCE POSSIBLE )
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/as_mlclci.h"
#include "asterfort/as_mlclor.h"
#include "asterfort/assert.h"
#include "asterfort/elraga.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/irnbsp.h"
!
    integer(kind=8) :: tygeos, nbpgm, nbpga, nloc, nutyma
    med_idt :: idfimd
    integer(kind=8) :: permu(nbpgm), codret, nbsp
    character(len=8) :: elrefa, fapg, nomtm
    character(len=64) :: locnam
!
    integer(kind=8) :: vali(2), jcopga, jwpga, ndim, nbfpg, iloc, dime, nbpgm2
    integer(kind=8) :: typgeo, nbpg, iret, nnoref, npgref, jrefco, jgscoo, jwg, jcorre
    integer(kind=8) :: ncorre, igau, idim, ad, ipgm, ipga, ada, im, ifm, nivinf
    character(len=8), parameter :: valk(3) = ['X', 'Y', 'Z']
    character(len=64) :: locnam2, nomasu
    integer(kind=8) :: edfuin
    parameter(edfuin=0)
    real(kind=8) :: xpgm, ypgm, zpgm, xpga, ypga, zpga, valr(2)
!
    call jemarq()
!
    call infniv(ifm, nivinf)
!   Voir irmpg1 pour la localisation des informations
    if (locnam(33:40) .eq. 'ASTER_SP') then
        call irnbsp(locnam(17:24), nbsp)
        nbpgm2 = nbpgm/nbsp

    else
        nbsp = 1
        nbpgm2 = nbpgm
    end if
!
!     DETERMINATION DES COORDONNES DES PG
!     DE L'ELEMENT DE REFERENCE ASTER
!     -------------------------------
    call wkvect('&&LRVCPG_COORD_PG_ASTER', 'V V R', 3*nbpga, jcopga)
    call wkvect('&&LRVCPG_POIDS_PG_ASTER', 'V V R', nbpga, jwpga)
    call elraga(elrefa, fapg, dime, nbfpg, zr(jcopga), &
                zr(jwpga))
!
!     NUMERO DU TYPE DE MAILLE DE L'ELEMENT DE REFERENCE : NUTYMA
    call jenonu(jexnom('&CATA.TM.NOMTM', nomtm), nutyma)
!
!     VERIFICATION SUR LE NOMBRE DE POINTS DE GAUSS
!     ---------------------------------------------
!     ATTENTION: POUR LES CHAMPS A SOUS-POINTS,
!     LE NBRE DE PG D'UN ELEMENT DE REF MED
!     PREND EN COMPTE LE NOMBRE DE SOUS-POINTS
!     CECI PEUT ETRE LA CAUSE DE L'EMISSION
!     DU MESSAGE CI-DESSOUS
    if (nbpgm2 .ne. nbpga) then
        vali(1) = nbpgm
        vali(2) = nbpga
        call utmess('A', 'MED_2', ni=2, vali=vali)
        codret = 4
        goto 999
    end if
!
!     DETERMINATION DU NOM DE LA LOCALISATION DES PG PRESENTE
!     DANS LE FICHIER MED ET CORRESPONDANT A L'ELEM DE REF ASTER.
!     ----------------------------------------------------------
!     -SI LA LOCALISATION EST ABSENTE : ON PREND EN COMPTE LA
!      LOCALISATION ASTER --> RISQUE DE RESULTATS FAUX
!     -SI LA LOCALISATION EST PRESENTE, ON COMPARE LES
!      COORDONNES DES PG ASTER/MED
    if (nivinf .gt. 1) then
        write (ifm, 101) nloc
    end if
    do iloc = 1, nloc
        call as_mlclci(idfimd, iloc, locnam2, tygeos, nbpg, &
                       ndim, nomasu, iret)
        if (tygeos .eq. typgeo) then
            if (locnam .eq. ' ' .or. locnam .eq. locnam2) goto 140
        end if
    end do
!     SI ON EST ICI, CELA SIGNIFIE QU'AUCUNE LOCALISATION
!     N'A ETE IDENTIFIEE POUR L'ELEMENT DE REFERENCE EN COURS
    if (nbpga .ne. 1 .or. nbpgm .ne. 1) then
        call utmess('A', 'MED_1', sk=elrefa)
    end if
    codret = 2
    goto 999
140 continue
!
!
    if (nivinf .gt. 1) then
        write (ifm, 102) locnam
    end if
!
!     DETERMINATION DES COORDONNES DES PG
!     DE L'ELEMENT DE REFERENCE MED
!     -------------------------------
    nnoref = (typgeo/100)*mod(typgeo, 100)
    npgref = (typgeo/100)*nbpgm
    call wkvect('&&LRVCPG_COORD_NO_MED', 'V V R', nnoref, jrefco)
    call wkvect('&&LRVCPG_COORD_PG_MED', 'V V R', npgref, jgscoo)
    call wkvect('&&LRVCPG_POIDS_PG_MED', 'V V R', nbpgm, jwg)
    call as_mlclor(idfimd, zr(jrefco), zr(jgscoo), zr(jwg), edfuin, &
                   locnam, iret)
    ASSERT(typgeo/100 .eq. dime)
!
!     COMPARAISON DES COORD DES PG ENTRE ASTER ET MED
!     -----------------------------------------------
!     NOMBRE DE PG NON APPARENTES : NCORRE
!     TABLEAU DE TRAVAIL ZI(JCORRE) DIMENSIONNE AU NBRE DE PG QUI VAUT:
!        -LE NUMERO DU PG LOCAL SI LA CORRESPONDANCE N'A PAS EU LIEU
!        -0 SINON
!
    call wkvect('&&LRVCPG_CORRESP_PG', 'V V I', nbpgm2, jcorre)
    ncorre = 0
    do igau = 1, nbpgm2
        zi(jcorre+igau-1) = 0
        do idim = 1, dime
            ad = dime*nbsp*(igau-1)+idim
            ada = dime*(igau-1)+idim
            if (nivinf .gt. 1) then
                write (ifm, 110) igau, idim, zr(jgscoo+ad-1), zr(jcopga+ &
                                                                 ada-1)
            end if
            if (abs(zr(jgscoo+ad-1)-zr(jcopga+ada-1)) .gt. 1.d-3) then
                ncorre = ncorre+1
                zi(jcorre+igau-1) = igau
                exit
            end if
        end do
    end do
!
!     SI LES PG ASTER/MED CORRESPONDENT : TOUT VA BIEN
    if (ncorre .eq. 0) then
        codret = 0
        goto 999
    else
!        .. SINON, ON RECHERCHE UNE EVENTUELLE PERMUTATION:
!        PERMU = LE TABLEAU DE PERMUTATIONS DIMENSIONNE
!        AU NBRE DE PG: PERMU(NUM_PG_MED)=NUM_PG_ASTER
        do ipgm = 1, nbpgm2
            permu(ipgm) = 0
            if (zi(jcorre+ipgm-1) .eq. 0) then
                permu(ipgm) = ipgm
            else
                ad = dime*nbsp*(igau-1)
                xpgm = zr(jgscoo+ad+1-1)
                ypgm = 0.d0
                zpgm = 0.d0
                if (dime .ge. 2) ypgm = zr(jgscoo+ad+2-1)
                if (dime .ge. 3) zpgm = zr(jgscoo+ad+3-1)
                do ipga = 1, nbpgm2
                    ada = dime*(ipga-1)
                    xpga = zr(jcopga+ada+1-1)
                    ypga = 0.d0
                    zpga = 0.d0
                    if (dime .ge. 2) ypga = zr(jcopga+ada+2-1)
                    if (dime .ge. 3) zpga = zr(jcopga+ada+3-1)
                    if (abs(xpgm-xpga) .lt. 1.d-3 .and. abs(ypgm-ypga) .lt. 1.d-3 .and. &
                        abs(zpgm-zpga) .lt. 1.d-3) then
                        permu(ipgm) = ipga
                        codret = 1
                        exit
                    end if
!                 SI ON EST ICI, CELA SIGNIFIE QUE L'UN DES PG MED
!                 N'A PAS PU ETRE IDENTIFIE A L'UN DES PG ASTER
!                 --> INCOMPATIBILITE DES PG, RISQUE DE RESULTATS FAUX
                    if (ipga .eq. nbpgm2) then
                        do im = 1, nbpgm2
                            call utmess('A+', 'MED_4', si=im)
                            do idim = 1, dime
                                valr(1) = zr(jgscoo+dime*nbsp*(im-1)+idim-1)
                                valr(2) = zr(jcopga+dime*(im-1)+idim-1)
                                call utmess('A+', 'MED_5', sk=valk(idim), nr=2, valr=valr)
                            end do
                        end do
                        call utmess('A', 'MED_3')
                        codret = 2
                        goto 999
                    end if
                end do
            end if
        end do
    end if
!
    if (codret .eq. 1) then
!        AFFICHAGE DES COORD DES PG MED/ASTER POUR
!        METTRE EN EVIDENCE LES PERMUTATIONS
        do im = 1, nbpgm2
            call utmess('A+', 'MED_4', si=im)
            do idim = 1, dime
                valr(1) = zr(jgscoo+dime*nbsp*(im-1)+idim-1)
                valr(2) = zr(jcopga+dime*(im-1)+idim-1)
                call utmess('A+', 'MED_5', sk=valk(idim), nr=2, valr=valr)
            end do
        end do
        call utmess('A', 'MED_6')
    end if
!
999 continue
!
    call jedetr('&&LRVCPG_COORD_PG_ASTER')
    call jedetr('&&LRVCPG_POIDS_PG_ASTER')
    call jedetr('&&LRVCPG_COORD_NO_MED')
    call jedetr('&&LRVCPG_COORD_PG_MED')
    call jedetr('&&LRVCPG_POIDS_PG_MED')
    call jedetr('&&LRVCPG_CORRESP_PG')
!
    call jedema()
!
101 format('  NOMBRE DE LOCALISATIONS LUES :', i4)
102 format('  LOCALISATION MED UTILISEE    :', a32)
110 format('  PT GAUSS', i4, ' DIM ', i1, ' COORD MED', 1pe12.5,&
     &                                  ' COORD ASTER', 1pe12.5)
end subroutine
