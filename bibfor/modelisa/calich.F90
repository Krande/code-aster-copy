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

subroutine calich(chargz, phenom)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/afrela.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/exisd.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=*), intent(in) :: chargz
    character(len=4), intent(in) :: phenom
!
!       CALICH -- TRAITEMENT DU MOT FACTEUR LIAISON_CHAMNO
!
!      TRAITEMENT DU MOT FACTEUR LIAISON_CHAMNO DE AFFE_CHAR_MECA
!      CE MOT FACTEUR PERMET DE DEFINIR UNE RELATION LINEAIRE ENTRE
!      LES DDLS DES NOEUDS D'UN MODELE DONT LES COEFFICIENTS SONT
!      LES VALEURS DES COMPOSANTES DU CHAM_NO DONNE APRES LE MOT CLE :
!      CHAM_NO.
!      LA VALEUR DU SECOND MEMBRE EST DONNEE APRES LE MOT CLE
!      COEF_IMPO (C'EST UN REEL).
!      ON NE PREND EN COMPTE QUE LES COEFFICIENTS NON NULS DU
!      CHAM_NO DANS LA RELATION LINEAIRE.
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
    character(len=19) :: lisrel, chamno, numeq
    complex(kind=8) :: betac
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    integer(kind=8) ::    iequa, ino, inocmp
    integer(kind=8) :: iocc, iret, k, nb, nbcmp, nbec, nbnoeu
    integer(kind=8) :: nbterm, nequa, nliai, nucmp
    real(kind=8) :: beta, vale, zero
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    real(kind=8), pointer :: vvale(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    motfac = 'LIAISON_CHAMNO'
!
    call getfac(motfac, nliai)
    if (nliai .eq. 0) goto 40
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
    lisrel = '&CALICH.RLLISTE'
!
! --- BOUCLE SUR LES OCCURENCES DU MOT-FACTEUR LIAISON_CHAMNO :
!     -------------------------------------------------------
    do iocc = 1, nliai
!
! ---   RECUPERATION DU CHAMNO
!       ----------------------
        call getvid(motfac, 'CHAM_NO', iocc=iocc, scal=chamno, nbret=nb)
        if (nb .eq. 0) then
            call utmess('F', 'CHARGES2_10')
        end if
!
! ---   VERIFICATION DE L'EXISTENCE DU CHAMNO
!       -------------------------------------
        call exisd("CHAM_NO", chamno, iret)
        if (iret .eq. 0) then
            call utmess('F', 'CHARGES2_11')
        end if
!
! ---   VERIFICATION DU TYPE DU CHAMP
!       -----------------------------
        call dismoi('TYPE_CHAMP', chamno, 'CHAM_NO', repk=tych)
        ASSERT(tych .eq. 'NOEU')
!
! ---   RECUPERATION DE LA VALEUR DU SECOND MEMBRE DE LA RELATION
! ---   LINEAIRE
!       --------
        call getvr8(motfac, 'COEF_IMPO', iocc=iocc, scal=beta, nbret=nb)
        ASSERT(nb .ne. 0)
!
! ---   RECUPERATION DE LA GRANDEUR ASSOCIEE AU CHAMNO :
!       ----------------------------------------------
        call dismoi('NOM_GD', chamno, 'CHAM_NO', repk=nomgd)
!
! ---   RECUPERATION DU NOMBRE DE MOTS SUR-LESQUELS SONT CODEES LES
! ---   LES INCONNUES ASSOCIEES A LA GRANDEUR DE NOM NOMGD
!       --------------------------------------------------
        call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nbec)
        ASSERT(nbec .le. 10)
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
! ---   DETERMINATION DU NOMBRE DE COMPOSANTES NON-NULLES DU CHAM_NO
!       ------------------------------------------------------------
        k = 0
        do i = 1, nequa
            if (vvale(i) .ne. zero) then
                k = k+1
            end if
        end do
!
        nbterm = k

        if (nbterm .eq. 0) then
            call utmess('F', 'CHARGES2_12')
        end if
!
! ---   CREATION DES TABLEAUX DE TRAVAIL NECESSAIRES A L'AFFECTATION
! ---   DE LA LISTE_RELA
!       ----------------
! ---     VECTEUR DU NOM DES NOEUDS
        AS_ALLOCATE(vk8=lisno, size=nbterm)
! ---     VECTEUR DU NOM DES DDLS
        AS_ALLOCATE(vk8=lisddl, size=nbterm)
! ---      VECTEUR DES COEFFICIENTS REELS
        AS_ALLOCATE(vr=coer, size=nbterm)
! ---     VECTEUR DES COEFFICIENTS COMPLEXES
        AS_ALLOCATE(vc=coec, size=nbterm)
! ---     VECTEUR DES DIRECTIONS DES DDLS A CONTRAINDRE
        AS_ALLOCATE(vr=direct, size=3*nbterm)
! ---     VECTEUR DES DIMENSIONS DE CES DIRECTIONS
        AS_ALLOCATE(vi=dime, size=nbterm)
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
                if (vale .ne. zero) then
                    k = k+1
                    nomcmp = zk8(inocmp+nucmp-1)
                    lisno(k) = nomnoe
                    lisddl(k) = nomcmp
                    coer(k) = vale
                end if
            end if
!
        end do
!
        nbterm = k
!
! ---   AFFECTATION DE LA RELATION A LA LISTE_RELA  :
!       ------------------------------------------
        call afrela(coer, coec, lisddl, lisno, dime, &
                    direct, nbterm, beta, betac, betaf, &
                    typcoe, typval, 0.d0, lisrel)
!
! ---   MENAGE :
!       ------
        AS_DEALLOCATE(vk8=lisno)
        AS_DEALLOCATE(vk8=lisddl)
        AS_DEALLOCATE(vr=coer)
        AS_DEALLOCATE(vc=coec)
        AS_DEALLOCATE(vr=direct)
        AS_DEALLOCATE(vi=dime)
!
    end do
!
! --- AFFECTATION DE LA LISTE_RELA A LA CHARGE :
!     ----------------------------------------
    if (phenom .eq. 'MECA') then
    end if
    call aflrch(lisrel, charge, 'LIN')
!
! --- MENAGE :
!     ------
    call jedetr(lisrel)
!
40  continue
!
    call jedema()
!.============================ FIN DE LA ROUTINE ======================
end subroutine
