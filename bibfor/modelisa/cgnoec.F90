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

subroutine cgnoec(mofaz, iocc, nomaz, lisnoz, nbno)
!.======================================================================
! aslint: disable=
    implicit none
!
!       CGNOEC -- TRAITEMENT DE L'OPTION ENV_CYLINDRE
!                 DU MOT FACTEUR CREA_GROUP_NO DE
!                 LA COMMANDE DEFI_GROUP
!
!      CETTE FONCTIONNALITE PERMET DE CREER UN GROUP_NO CONSTITUE
!      DE TOUS LES NOEUDS APPARTENANT A L'ENVELOPPE D'UN CYLINDRE
!      DEFINI PAR L'UTILISATEUR.
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_NO'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISNOZ        - JXVAR - K24  - : NOM DE LA LISTE DE NOEUDS
!                                   APPARTENANT A L'ENVELOPPE
!                                   DU CYLINDRE.
!  NBNO          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!.========================= DEBUT DES DECLARATIONS ====================
!
! -----  ARGUMENTS
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utcono.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: mofaz, nomaz, lisnoz
!
! --------- VARIABLES LOCALES ---------------------------
    character(len=8) :: noma, k8bid
    character(len=16) :: motfac, mocle(3)
    character(len=24) :: lisnoe
!
    real(kind=8) :: x0(3), x(3), xx0(3), axe(3), angle(2)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!-----------------------------------------------------------------------
    integer(kind=8) ::  idlino, ino, iocc, iret, nangle
    integer(kind=8) :: nb, nbno, nbnoe, ndim, nprec, nrayon, nv
    integer(kind=8) :: nvect
    real(kind=8) :: ang, d2, dist, eps, prec, psca
    real(kind=8) :: rayon, un, xnorm, xnorm2, xnoxx0
    real(kind=8) :: xnoxx2, zero
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    motfac = mofaz
    noma = nomaz
    lisnoe = lisnoz
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
    rayon = zero
!
    eps = 100.0d0*r8prem()
!
    nbno = 0
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
        call utmess('F', 'MODELISA3_84')
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
            call utmess('F', 'MODELISA3_85')
        else
            nvect = -nvect
            if (nvect .ne. 3) then
                call utmess('F', 'MODELISA3_86')
            else
                call getvr8(motfac, 'VECT_NORMALE', iocc=iocc, nbval=nvect, vect=axe, &
                            nbret=nv)
            end if
        end if
    else
        nangle = -nangle
        if (nangle .ne. 2) then
            call utmess('F', 'MODELISA3_87')
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
! --- RECUPERATION DE LA DEMI-EPAISSEUR DE L'ENVELOPPE :
!     ------------------------------------------------
    call getvr8(motfac, 'PRECISION', iocc=iocc, nbval=0, nbret=nprec)
    if (nprec .eq. 0) then
        call utmess('F', 'MODELISA3_88')
    else
        call getvr8(motfac, 'PRECISION', iocc=iocc, scal=prec, nbret=nb)
        if (prec .le. zero) then
            call utmess('F', 'MODELISA3_89')
        end if
    end if
!
! --- RECUPERATION DU NOMBRE DE NOEUDS DU MAILLAGE :
!     ---------------------------------------------
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoe)
!
! --- ALLOCATION DU VECTEUR DES NOMS DES NOEUDS  APPARTENANT
! --- A L'ENVELOPPE DU CYLINDRE :
!     -------------------------
    call wkvect(lisnoe, 'V V I', nbnoe, idlino)
!
! --- PARCOURS DES NOEUDS DU MAILLAGE :
!     --------------------------------
    nbno = 0
    do ino = 1, nbnoe
!
! ---     COORDONNEES DU NOEUD :
!         --------------------
        x(1) = vale(3*(ino-1)+1)
        x(2) = vale(3*(ino-1)+2)
        x(3) = vale(3*(ino-1)+3)
!
        xx0(1) = x(1)-x0(1)
        xx0(2) = x(2)-x0(2)
        xx0(3) = x(3)-x0(3)
!
        xnoxx2 = xx0(1)*xx0(1)+xx0(2)*xx0(2)+xx0(3)*xx0(3)
!
        if (xnoxx2 .ne. zero) then
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
                 x(1)-x0(1))*(x(1)-x0(1))+(x(2)-x0(2))*(x(2)-x0(2))+(x(3)-x0(3))*(x(3)-x0(3)&
                 &))*sin(ang &
                 )*sin(ang &
                 )
!
! ---         SI LE NOEUD COURANT APPARTIENT A L'ENVELOPPE DU
! ---         CYLINDRE, ON L'AFFECTE A LA LISTE DE NOEUDS QUI
! ---         SERA AFFECTEE AU GROUP_NO :
!             -------------------------
            dist = sqrt(d2)
            if (abs(dist-rayon) .le. prec) then
                nbno = nbno+1
                zi(idlino+nbno-1) = ino
            end if
!
        else
!              -- CAS DU NOEUD CONFONDU AVEC LE POINT DEFINISSANT L'AXE:
            if (abs(rayon) .le. prec) then
                nbno = nbno+1
                zi(idlino+nbno-1) = ino
            end if
        end if
!
    end do
!
    call jedema()
!.============================ FIN DE LA ROUTINE ======================
end subroutine
