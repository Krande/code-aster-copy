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

subroutine cgmasp(mofaz, iocc, nomaz, lismaz, nbma)
!.======================================================================
    implicit none
!
!       CGMASP -- TRAITEMENT DE L'OPTION SPHERE
!                 DU MOT FACTEUR CREA_GROUP_MA DE
!                 LA COMMANDE DEFI_GROUP
!
!      CETTE FONCTIONNALITE PERMET DE CREER UN GROUP_MA CONSTITUE
!      DE TOUTES LES MAILLES DONT UN NOEUD AU MOINS APPARTIENT
!      A UNE SPHERE DE RAYON R ET DE CENTRE P0 (X0,Y0,Z0).
!      LE RAYON ET LE POINT P0 SONT DES DONNEES UTILISATEUR.
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_MA'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISMAZ        - JXVAR - K24  - : NOM DE LA LISTE DE MAILLES DONT UN
!                                   NOEUD AU MOINS APPARTIENT A LA
!                                   SPHERE DEFINIE PAR L'UTILISATEUR
!  NBMA          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!.========================= DEBUT DES DECLARATIONS ====================
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utcono.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
!
! -----  ARGUMENTS
    character(len=*) :: mofaz, nomaz, lismaz
!
! --------- VARIABLES LOCALES ---------------------------
    integer(kind=8) :: nbnod, nbno
    character(len=8) :: noma, k8bid, nomail
    character(len=16) :: motfac, mocle(3)
    character(len=16) :: selec
    character(len=24) :: lismai
!
    real(kind=8) :: x0(3), x(3)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, idlima, idnoeu, ima, ino
    integer(kind=8) :: iocc, iret, nb, nbma, nbmai, ndim, nrayon
    integer(kind=8) :: numnoe
    real(kind=8) :: d2, rayon, zero
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    motfac = mofaz
    noma = nomaz
    lismai = lismaz
!
! --- RECUPERATION DU TYPE DE VERIFICATION A APPLIQUER :
!     --------------------------------------------------
    call getvtx(motfac, 'CRIT_NOEUD', iocc=iocc, scal=selec, nbret=ibid)
!
    zero = 0.0d0
!
    x0(1) = zero
    x0(2) = zero
    x0(3) = zero
!
    x(1) = zero
    x(2) = zero
    x(3) = zero
!
    rayon = zero
!
    nbma = 0
!
! --- RECUPERATION DE LA DIMENSION DU MAILLAGE :
!     ----------------------------------------
    call dismoi('Z_CST', noma, 'MAILLAGE', repk=k8bid)
    if (k8bid(1:3) .eq. 'OUI') then
        ndim = 2
    else
        ndim = 3
    end if
!
! --- RECUPERATION DES COORDONNES DES NOEUDS DU MAILLAGE :
!     --------------------------------------------------
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! --- RECUPERATION DU CENTRE DE LA SPHERE (OU DU CERCLE) :
!     --------------------------------------------------
    mocle(1) = 'POINT'
    mocle(2) = 'NOEUD_CENTRE'
    mocle(3) = 'GROUP_NO_CENTRE'
    call utcono(motfac, mocle, iocc, noma, ndim, &
                x0, iret)
!
! --- RECUPERATION DU RAYON DE LA SPHERE :
!     ----------------------------------
    call getvr8(motfac, 'RAYON', iocc=iocc, nbval=0, nbret=nrayon)
    if (nrayon .eq. 0) then
        call utmess('F', 'MODELISA3_82')
    else
        call getvr8(motfac, 'RAYON', iocc=iocc, scal=rayon, nbret=nb)
        if (rayon .le. zero) then
            call utmess('F', 'MODELISA3_83')
        end if
    end if
!
! --- RECUPERATION DU NOMBRE DE MAILLES DU MAILLAGE :
!     ---------------------------------------------
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmai)
!
! --- ALLOCATION DU VECTEUR DES NOMS DES MAILLES  APPARTENANT
! --- A LA SPHERE :
!     -----------
    call wkvect(lismai, 'V V I', nbmai, idlima)
!
! --- PARCOURS DES MAILLES DU MAILLAGE :
!     --------------------------------
    do ima = 1, nbmai
!
! ---     RECUPERATION DU NOM DE LA MAILLE à partrir du numero d'ordre:
!         --------------------------------
        nomail = int_to_char8(ima)
!
! ---     RECUPERATION DES CONNECTIVITES DE LA MAILLE :
!         -------------------------------------------
        ibid = char8_to_int(nomail)

        call jeveuo(jexnum(noma//'.CONNEX', ibid), 'L', idnoeu)
!
! ---     RECUPERATION DU NOMBRE DE CONNECTIVITES DE LA MAILLE :
!         ----------------------------------------------------
        ibid = char8_to_int(nomail)
        call jelira(jexnum(noma//'.CONNEX', ibid), 'LONMAX', nbno)
!
! ---      COMPTE NOMBRE DES NOEUDS D'UN MAILLE DANS LE SPHERE :
!          ----------------------------------------------------
        nbnod = 0
!
! ---     BOUCLE SUR LES CONNECTIVITES DE LA MAILLE :
!         -----------------------------------------
        do ino = 1, nbno
!
! ---        NUMERO DU NOEUD :
!            ---------------
            numnoe = zi(idnoeu+ino-1)
!
! ---        COORDONNEES DU NOEUD :
!            --------------------
            x(1) = vale(3*(numnoe-1)+1)
            x(2) = vale(3*(numnoe-1)+2)
            if (ndim .eq. 3) then
                x(3) = vale(3*(numnoe-1)+3)
            end if
!
! ---        DISTANCE DU NOEUD COURANT AU CENTRE DE LA SPHERE :
!            ------------------------------------------------
            d2 = (x(1)-x0(1))*(x(1)-x0(1))+(x(2)-x0(2))*(x(2)-x0(2))+(x(3)-x0(3))*(x(3)-x0(3) &
                                                                                   )
!
! ---      SI LE MOT CLE SIMPLE CRIT_NOEUD EST EGAL A AU MOINS UN NOEUD
!          -------------------------------------------------------------
            if (selec .eq. 'AU_MOINS_UN') then
!
! ---             SI LE NOEUD COURANT EST DANS LA SPHERE, ON AFFECTE
! ---             LA MAILLE COURANTE A LA LISTE DE MAILLES QUI SERA
! ---             AFFECTEE AU GROUP_MA :
!                 --------------------
                if (d2 .le. rayon*rayon) then
                    nbma = nbma+1
                    zi(idlima+nbma-1) = ima
                    nomail = int_to_char8(ima)
                    goto 10
                end if
! ---            SI LE MOT CLE SIMPLE CRIT_NOEUD EST EGAL A TOUT OU
! ---            MAJORITE , COMPTER LE NOMBRE DES NOEUDS D'UNE MAILLE
! ---            DANS LE SPHERE :
!                ----------------------------------------------------
            else if ((selec .eq. 'TOUS') .or. (selec .eq. 'MAJORITE')) then
                if (d2 .le. rayon*rayon) then
                    nbnod = nbnod+1
                end if
            end if
!
        end do
!
        if (selec .eq. 'TOUS') then
            if (nbnod .eq. nbno) then
                nbma = nbma+1
                zi(idlima+nbma-1) = ima
                nomail = int_to_char8(ima)
                goto 10
            end if
        end if
        if (selec .eq. 'MAJORITE') then
            if (nbnod .ge. (nbno+1)/2) then
                nbma = nbma+1
                zi(idlima+nbma-1) = ima
                nomail = int_to_char8(ima)
                goto 10
            end if
        end if
!
10      continue
    end do
!
    call jedema()
!.============================ FIN DE LA ROUTINE ======================
end subroutine
