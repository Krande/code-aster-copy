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

subroutine cgmacy(mofaz, iocc, nomaz, lismaz, nbma)
!.======================================================================
! aslint: disable=
    implicit none
!
!       CGMACY -- TRAITEMENT DE L'OPTION CYLINDRE
!                 DU MOT FACTEUR CREA_GROUP_MA DE
!                 LA COMMANDE DEFI_GROUP
!
!      CETTE FONCTIONNALITE PERMET DE CREER UN GROUP_MA CONSTITUE
!      DE TOUTES LES MAILLES DONT UN NOEUD AU MOINS APPARTIENT
!      A CYLINDRE DONT L'AXE ET LE RAYON SONT DEFINIS PAR
!      L'UTILISATEUR.
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_MA'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISMAZ        - JXVAR - K24  - : NOM DE LA LISTE DE MAILLES DONT UN
!                                   NOEUD AU MOINS APPARTIENT AU
!                                   CYLINDRE DEFINI PAR L'UTILISATEUR
!  NBMA          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!.========================= DEBUT DES DECLARATIONS ====================
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
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
    character(len=8) :: noma, k8bid, nomail
    character(len=16) :: motfac, mocle(3)
    character(len=24) :: lismai
    character(len=16) :: selec
!
!
    real(kind=8) :: x0(3), x(3), xx0(3), axe(3), angle(2)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, idlima, idnoeu, ima, ino
    integer(kind=8) :: iocc, iret, nangle, nb, nbma, nbmai, nbno
    integer(kind=8) :: nbnod, ndim, nrayon, numnoe, nv, nvect
    real(kind=8) :: ang, d2, eps, psca
    real(kind=8) :: rayon, un, xnorm, xnorm2, xnoxx0, xnoxx2, zero
    real(kind=8), pointer :: vale(:) => null()
!
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
    un = 1.0d0
!
    x0(1) = zero
    x0(2) = zero
    x0(3) = zero
!
    x(1) = zero
    x(2) = zero
    x(3) = zero
!
    xx0(1) = zero
    xx0(2) = zero
    xx0(3) = zero
!
    axe(1) = zero
    axe(2) = zero
    axe(3) = zero
!
    rayon = zero
!
    eps = 100.0d0*r8prem()
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
    if (ndim .ne. 3) then
        call utmess('F', 'MODELISA3_73')
    end if
!
! --- RECUPERATION DES COORDONNES DES NOEUDS DU MAILLAGE :
!     --------------------------------------------------
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! --- RECUPERATION DU POINT SITUE SUR L'AXE DU CYLINDRE :
!     -------------------------------------------------
    mocle(1) = 'POINT'
    mocle(2) = 'NOEUD_CENTRE'
    mocle(3) = 'GROUP_NO_CENTRE'
    call utcono(motfac, mocle, iocc, noma, ndim, &
                x0, iret)
!
! --- RECUPERATION DU RAYON DU CYLINDRE :
!     ---------------------------------
    call getvr8(motfac, 'RAYON', iocc=iocc, nbval=0, nbret=nrayon)
    if (nrayon .eq. 0) then
        call utmess('F', 'MODELISA3_74')
    else
        call getvr8(motfac, 'RAYON', iocc=iocc, scal=rayon, nbret=nb)
        if (rayon .le. zero) then
            call utmess('F', 'MODELISA3_75')
        end if
    end if
!
! --- RECUPERATION DE LA DIRECTION DEFINISSANT L'AXE DU CYLINDRE :
!     ----------------------------------------------------------
    call getvr8(motfac, 'ANGL_NAUT', iocc=iocc, nbval=0, nbret=nangle)
    if (nangle .eq. 0) then
        call getvr8(motfac, 'VECT_NORMALE', iocc=iocc, nbval=0, nbret=nvect)
        if (nvect .eq. 0) then
            call utmess('F', 'MODELISA3_76')
        else
            nvect = -nvect
            if (nvect .ne. 3) then
                call utmess('F', 'MODELISA3_77')
            else
                call getvr8(motfac, 'VECT_NORMALE', iocc=iocc, nbval=nvect, vect=axe, &
                            nbret=nv)
            end if
        end if
    else
        nangle = -nangle
        if (nangle .ne. 2) then
            call utmess('F', 'MODELISA3_78')
        end if
        call getvr8(motfac, 'ANGL_NAUT', iocc=iocc, nbval=nangle, vect=angle, &
                    nbret=nv)
!
        angle(1) = angle(1)*r8dgrd()
        angle(2) = angle(2)*r8dgrd()
!
        axe(1) = cos(angle(1))*cos(angle(2))
        axe(2) = sin(angle(1))*cos(angle(2))
        axe(3) = -sin(angle(2))
    end if
!
    xnorm2 = axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3)
!
    if (xnorm2 .eq. zero) then
        call utmess('F', 'MODELISA3_79')
    end if
!
    xnorm = sqrt(xnorm2)
!
    axe(1) = axe(1)/xnorm
    axe(2) = axe(2)/xnorm
    axe(3) = axe(3)/xnorm
!
! --- RECUPERATION DU NOMBRE DE MAILLES DU MAILLAGE :
!     ---------------------------------------------
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmai)
!
! --- ALLOCATION DU VECTEUR DES NOMS DES MAILLES APPARTENANT AU
! --- CYLINDRE :
!     --------
    call wkvect(lismai, 'V V I', nbmai, idlima)
!
! --- PARCOURS DES MAILLES DU MAILLAGE :
!     --------------------------------
    do ima = 1, nbmai
!
! ---     RECUPERATION DU NOM DE LA MAILLE :
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
! ---      COMPTE NOMBRE DES NOEUDS D'UN MAILLE DANS LE CYLINDRE :
!          ------------------------------------------------------
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
            x(3) = vale(3*(numnoe-1)+3)
!
            xx0(1) = x(1)-x0(1)
            xx0(2) = x(2)-x0(2)
            xx0(3) = x(3)-x0(3)
!
            xnoxx2 = xx0(1)*xx0(1)+xx0(2)*xx0(2)+xx0(3)*xx0(3)
!
! ---       SI LE NOEUD COURANT EST CONFONDU AVEC LE NOEUD DE
! ---       L'AXE DU CYLINDRE, ON AFFECTE  LA MAILLE COURANTE A LA
! ---       LISTE DE MAILLES QUI SERA AFFECTEE AU GROUP_MA :
!           ----------------------------------------------
            if (xnoxx2 .eq. zero) then
                nbma = nbma+1
                zi(idlima+nbma-1) = ima
                goto 10
            else
!
                xnoxx0 = sqrt(xnoxx2)
!
                xx0(1) = xx0(1)/xnoxx0
                xx0(2) = xx0(2)/xnoxx0
                xx0(3) = xx0(3)/xnoxx0
!
! ---         CALCUL DE L'ANGLE FORME PAR L'AXE DU CYLINDRE
! ---         AVEC LE VECTEUR POSITION COURANT XX0 :
!             ------------------------------------
                psca = abs(xx0(1)*axe(1)+xx0(2)*axe(2)+xx0(3)*axe(3))
                if (psca .gt. un) then
                    psca = psca-eps
                end if
                ang = acos(psca)
!
! ---         CALCUL DE LA DISTANCE DU NOEUD COURANT A L'AXE
! ---         DU CYLINDRE :
!             -----------
                d2 = ( &
                     ( &
                     x(1)-x0(1))*(x(1)-x0(1))+(x(2)-x0(2))*(x(2)-x0(2))+(x(3)-x0(3))*(x(3)-x&
                     &0(3)))*sin(ang &
                     )*sin(ang &
                     )
!
! ---      SI LE MOT CLE SIMPLE CRIT_NOEUD EST EGAL A AU MOINS UN NOEUD
!          -------------------------------------------------------------
!
                if (selec .eq. 'AU_MOINS_UN') then
!
! ---         SI LE NOEUD COURANT EST DANS LE CYLINDRE, ON AFFECTE
! ---         LA MAILLE COURANTE A LA LISTE DE MAILLES QUI SERA
! ---         AFFECTEE AU GROUP_MA :
!             --------------------
                    if (d2 .le. rayon*rayon) then
                        nbma = nbma+1
                        zi(idlima+nbma-1) = ima
                        goto 10
                    end if
!
! ---         SI LE MOT CLE SIMPLE CRIT_NOEUD EST EGAL A TOUT OU
! ---         MAJORITE, COMPTER LE NOMBRE DES NOEUDS D'UNE MAILLE
! ---         DANS LE CYLINDRE :
!             ------------------------------------------------
                else if ((selec .eq. 'TOUS') .or. (selec .eq. 'MAJORITE')) &
                    then
!
                    if (d2 .le. rayon*rayon) then
                        nbnod = nbnod+1
                    end if
!
                end if
!
            end if
!
        end do
!
        if (selec .eq. 'TOUS') then
            if (nbnod .eq. nbno) then
                nbma = nbma+1
                zi(idlima+nbma-1) = ima
                goto 10
            end if
        end if
!
        if (selec .eq. 'MAJORITE') then
            if (nbnod .ge. (nbno+1)/2) then
                nbma = nbma+1
                zi(idlima+nbma-1) = ima
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
