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

subroutine caimch(chargz)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
    character(len=*) :: chargz
!
!       CAIMCH -- TRAITEMENT DU MOT FACTEUR CHAMNO_IMPO
!
!      TRAITEMENT DU MOT FACTEUR CHAMNO_IMPO DE AFFE_CHAR_MECA
!      CE MOT FACTEUR PERMET D'IMPOSER SUR DES DDL DES NOEUDS
!      D'UN MODELELES VALEURS DES COMPOSANTES DU CHAM_NO DONNE
!      APRES LE MOT CLE : CHAM_NO.
!
! -------------------------------------------------------
!  CHARGE        - IN    - K8   - : NOM DE LA SD CHARGE
!                - JXVAR -      -   LA  CHARGE EST ENRICHIE
!                                   DE LA RELATION LINEAIRE DECRITE
!                                   CI-DESSUS.
! -------------------------------------------------------
!
    character(len=4) :: tych, typval, typcoe
    character(len=8) :: noma, nomcmp, nomnoe, betaf
    character(len=8) :: charge, nomgd
    character(len=16) :: motfac
    character(len=19) :: lisrel, numeq, chamno
    real(kind=8) :: beta, coef_impo
    complex(kind=8) :: betac
    integer(kind=8) :: idcoec, idcoer, idddl, idimen, idirec
    integer(kind=8) :: idnoeu, iequa, ino, inocmp, iocc
    integer(kind=8) :: iret, k, nb, nbcmp, nbec, nbnoeu, nbterm
    integer(kind=8) :: nequa, nliai, nucmp
    real(kind=8) :: vale, zero
    integer(kind=8), pointer :: deeq(:) => null()
    real(kind=8), pointer :: vvale(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    motfac = 'CHAMNO_IMPO'
!
    call getfac(motfac, nliai)
    if (nliai .eq. 0) goto 30
!
! --- INITIALISATIONS :
!     ---------------
    zero = 0.0d0
!
! --- BETA, BETAC ET BETAF SONT LES VALEURS DU SECOND MEMBRE DE LA
! --- RELATION LINEAIRE SUIVANT QUE C'EST UN REEL, UN COMPLEXE OU
! --- UNE FONCTION, DANS NOTRE CAS C'EST UN REEL
!
    beta = zero
    betac = (0.0d0, 0.0d0)
    betaf = '&FOZERO'
!
    charge = chargz
!
! --- TYPE DES VALEURS AU SECOND MEMBRE DE LA RELATION
!
    typval = 'REEL'
!
! --- TYPE DES VALEURS DES COEFFICIENTS
!
    typcoe = 'REEL'
!
! --- NOM DE LA LISTE_RELA
!
    lisrel = '&CAIMCH.RLLISTE'
!
! --- BOUCLE SUR LES OCCURENCES DU MOT-FACTEUR CHAMNO_IMPO :
!     -------------------------------------------------------
    do iocc = 1, nliai
!
! ---   RECUPERATION DU CHAMNO
!       ----------------------
        call getvid(motfac, 'CHAM_NO', iocc=iocc, scal=chamno, nbret=nb)
        if (nb .eq. 0) then
            call utmess('F', 'MODELISA2_83')
        end if
!
! ---   VERIFICATION DE L'EXISTENCE DU CHAMNO
!       -------------------------------------
        call exisd("CHAM_NO", chamno, iret)
        if (iret .eq. 0) then
            call utmess('F', 'MODELISA2_84')
        end if
!
! ---   VERIFICATION DU TYPE DU CHAMP
!       -----------------------------
        call dismoi('TYPE_CHAMP', chamno, 'CHAM_NO', repk=tych)
!
        if (tych .ne. 'NOEU') then
            call utmess('F', 'MODELISA2_85')
        end if
!
! ---   RECUPERATION DE LA VALEUR DU SECOND MEMBRE DE LA RELATION
! ---   LINEAIRE
!       --------
        call getvr8(motfac, 'COEF_IMPO', iocc=iocc, scal=coef_impo, nbret=nb)
        if (nb .eq. 0) then
            call utmess('F', 'MODELISA2_86')
        end if
!
! ---   RECUPERATION DE LA GRANDEUR ASSOCIEE AU CHAMNO :
!       ----------------------------------------------
        call dismoi('NOM_GD', chamno, 'CHAM_NO', repk=nomgd)
!
! ---   RECUPERATION DU NOMBRE DE MOTS SUR-LESQUELS SONT CODEES LES
! ---   LES INCONNUES ASSOCIEES A LA GRANDEUR DE NOM NOMGD
!       --------------------------------------------------
        call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nbec)
        if (nbec .gt. 11) then
            call utmess('F', 'MODELISA2_87', sk=nomgd)
        end if
!
! ---   RECUPERATION DU MAILLAGE ASSOCIE AU CHAM_NO
!       -------------------------------------------
        call dismoi('NOM_MAILLA', chamno, 'CHAM_NO', repk=noma)
!
! ---   RECUPERATION DU NOMBRE DE NOEUDS DU MAILLAGE
!       --------------------------------------------
        call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoeu)
!
! ---   RECUPERATION DU NOMBRE DE TERMES DU CHAM_NO
!       -------------------------------------------
        call dismoi('NB_EQUA', chamno, 'CHAM_NO', repi=nequa)
!
! ---   RECUPERATION DU NUME_EQUA DU CHAM_NO
!       ------------------------------------
        call dismoi('NUME_EQUA', chamno, 'CHAM_NO', repk=numeq)
!
! ---   RECUPERATION DU NOMBRE DE COMPOSANTES ASSOCIEES A LA LA GRANDEUR
!       ----------------------------------------------------------------
        call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', nbcmp)
!
! ---   RECUPERATION DU NOM DES COMPOSANTES ASSOCIEES A LA LA GRANDEUR
!       --------------------------------------------------------------
        call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', inocmp)
!
! ---   RECUPERATION DU .VALE DU CHAM_NO
!       --------------------------------
        call jeveuo(chamno//'.VALE', 'E', vr=vvale)
!
! ---   RECUPERATION DU .DEEQ DU NUME_EQUA
!       ----------------------------------
        call jeveuo(numeq//'.DEEQ', 'L', vi=deeq)
!
        nbterm = nequa
!
! ---   CREATION DES TABLEAUX DE TRAVAIL NECESSAIRES A L'AFFECTATION
! ---   DE LA LISTE_RELA
!       ----------------
! ---     VECTEUR DU NOM DES NOEUDS
        call wkvect('&&CAIMCH.LISNO', 'V V K8', nbterm, idnoeu)
! ---     VECTEUR DU NOM DES DDLS
        call wkvect('&&CAIMCH.LISDDL', 'V V K8', nbterm, idddl)
! ---      VECTEUR DES COEFFICIENTS REELS
        call wkvect('&&CAIMCH.COER', 'V V R', nbterm, idcoer)
! ---     VECTEUR DES COEFFICIENTS COMPLEXES
        call wkvect('&&CAIMCH.COEC', 'V V C', nbterm, idcoec)
! ---     VECTEUR DES DIRECTIONS DES DDLS A CONTRAINDRE
        call wkvect('&&CAIMCH.DIRECT', 'V V R', 3*nbterm, idirec)
! ---     VECTEUR DES DIMENSIONS DE CES DIRECTIONS
        call wkvect('&&CAIMCH.DIME', 'V V I', nbterm, idimen)
!
! ---   AFFECTATION DES TABLEAUX DE TRAVAIL :
!       -----------------------------------
        k = 0
!
! ---   BOUCLE SUR LES TERMES DU CHAM_NO
!
        do iequa = 1, nequa
!
! ---     INO  : NUMERO DU NOEUD INO CORRESPONDANT AU DDL IEQUA
!
            ino = deeq(1+2*(iequa-1)+1-1)
!
! ---     NUCMP  : NUMERO DE COMPOSANTE CORRESPONDANTE AU DDL IEQUA
!
            nucmp = deeq(1+2*(iequa-1)+2-1)
!
! ---     ON NE PREND PAS EN COMPTE LES MULTIPLICATEURS DE LAGRANGE
! ---     (CAS OU NUCMP < 0)
!
            if (nucmp .gt. 0) then
!
! ---       RECUPERATION DU NOM DU NOEUD INO
!
                nomnoe = int_to_char8(ino)
!
                vale = vvale(iequa)
!
                k = k+1
                nomcmp = zk8(inocmp+nucmp-1)
                zk8(idnoeu+k-1) = nomnoe
                zk8(idddl+k-1) = nomcmp
                zr(idcoer+k-1) = 1.
                beta = coef_impo*vale
!
! ---       AFFECTATION DE LA RELATION A LA LISTE_RELA  :
!
                call afrela(zr(idcoer+k-1), zc(idcoec+k-1), zk8(idddl+k-1), zk8(idnoeu+k-1), &
                            zi(idimen+k-1), [0.d0], 1, beta, betac, &
                            betaf, typcoe, typval, 0.d0, lisrel)
            end if
!
        end do
!
        nbterm = k
!
! ---   MENAGE :
!       ------
        call jedetr('&&CAIMCH.LISNO')
        call jedetr('&&CAIMCH.LISDDL')
        call jedetr('&&CAIMCH.COER')
        call jedetr('&&CAIMCH.COEC')
        call jedetr('&&CAIMCH.DIRECT')
        call jedetr('&&CAIMCH.DIME')
!
    end do
!
! --- AFFECTATION DE LA LISTE_RELA A LA CHARGE :
!     ----------------------------------------
    call aflrch(lisrel, charge, 'LIN')
!
! --- MENAGE :
!     ------
    call jedetr(lisrel)
!
30  continue
!
    call jedema()
!
end subroutine
