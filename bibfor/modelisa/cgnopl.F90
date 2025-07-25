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

subroutine cgnopl(mofaz, iocc, nomaz, lisnoz, nbno)
!.======================================================================
    implicit none
!
!       CGNOPL -- TRAITEMENT DE L'OPTION PLAN
!                 DU MOT FACTEUR CREA_GROUP_NO DE
!                 LA COMMANDE DEFI_GROUP
!
!      CETTE FONCTIONNALITE PERMET DE CREER UN GROUP_NO CONSTITUE
!      DE TOUS LES NOEUDS APPARTENANT A UNE DROITE EN 2D
!      OU UN PLAN EN 3D DEFINIS PAR L'UTILISATEUR.
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_NO'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISNOZ        - JXVAR - K24  - : NOM DE LA LISTE DE NOEUDS
!                                   APPARTENANT A LA DROITE (EN 2D)
!                                   OU AU PLAN (EN 3D) DONNES PAR
!                                   L'UTILISATEUR
!  NBNO          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!.========================= DEBUT DES DECLARATIONS ====================
!
! -----  ARGUMENTS
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/cgnop0.h"
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
    real(kind=8) :: x0(3), vecnor(3), angle(2)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!-----------------------------------------------------------------------
    integer(kind=8) ::  idlino, iocc, iret, nangle
    integer(kind=8) :: nb, nbno, nbnoe, ndim, ndim1, nprec, nv
    integer(kind=8) :: nvect
    real(kind=8) :: prec, xnorm, xnorm2, zero
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
!
    x0(1) = zero
    x0(2) = zero
    x0(3) = zero
!
    vecnor(1) = zero
    vecnor(2) = zero
    vecnor(3) = zero
!
! --- RECUPERATION DE LA DIMENSION DU MAILLAGE :
!     ----------------------------------------
    ndim = 3
    call dismoi('Z_CST', noma, 'MAILLAGE', repk=k8bid)
    if (k8bid .eq. 'OUI') ndim = 2
!
! --- RECUPERATION DES COORDONNES DES NOEUDS DU MAILLAGE :
!     --------------------------------------------------
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! --- RECUPERATION DU POINT SITUE SUR LE PLAN OU LA DROITE :
!     ----------------------------------------------------
    mocle(1) = 'POINT'
    mocle(2) = 'NOEUD_CENTRE'
    mocle(3) = 'GROUP_NO_CENTRE'
    call utcono(motfac, mocle, iocc, noma, ndim, &
                x0, iret)
!
! --- RECUPERATION DE LA DIRECTION PERPENDICULAIRE AU PLAN MILIEU
! --- DE LA BANDE :
!     -----------
    call getvr8(motfac, 'ANGL_NAUT', iocc=iocc, nbval=0, nbret=nangle)
    if (nangle .eq. 0) then
        call getvr8(motfac, 'VECT_NORMALE', iocc=iocc, nbval=0, nbret=nvect)
        if (nvect .eq. 0) then
            call utmess('F', 'MODELISA3_93')
        else
            nvect = -nvect
            if (ndim .eq. 3 .and. nvect .ne. 3) then
                call utmess('F', 'MODELISA3_94')
            else if (ndim .eq. 2 .and. nvect .ne. 2) then
                call utmess('F', 'MODELISA3_95')
            else
                call getvr8(motfac, 'VECT_NORMALE', iocc=iocc, nbval=nvect, vect=vecnor, &
                            nbret=nv)
            end if
        end if
    else
        nangle = -nangle
        ndim1 = ndim-1
        nangle = min(nangle, ndim1)
        call getvr8(motfac, 'ANGL_NAUT', iocc=iocc, nbval=nangle, vect=angle, &
                    nbret=nv)
!
        if (ndim .eq. 2) then
            angle(1) = angle(1)*r8dgrd()
!
            vecnor(1) = cos(angle(1))
            vecnor(2) = sin(angle(1))
            vecnor(3) = zero
        else if (ndim .eq. 3) then
            angle(1) = angle(1)*r8dgrd()
            angle(2) = angle(2)*r8dgrd()
!
            vecnor(1) = cos(angle(1))*cos(angle(2))
            vecnor(2) = sin(angle(1))*cos(angle(2))
            vecnor(3) = -sin(angle(2))
        end if
    end if
!
    xnorm2 = vecnor(1)*vecnor(1)+vecnor(2)*vecnor(2)+vecnor(3)*vecnor(3)
!
    if (xnorm2 .eq. zero) then
        call utmess('F', 'MODELISA3_96')
    end if
!
    xnorm = sqrt(xnorm2)
!
    vecnor(1) = vecnor(1)/xnorm
    vecnor(2) = vecnor(2)/xnorm
    vecnor(3) = vecnor(3)/xnorm
!
! --- RECUPERATION DE LA TOLERANCE :
!     ----------------------------
    call getvr8(motfac, 'PRECISION', iocc=iocc, nbval=0, nbret=nprec)
    if (nprec .eq. 0) then
        call utmess('F', 'MODELISA3_97')
    else
        call getvr8(motfac, 'PRECISION', iocc=iocc, scal=prec, nbret=nb)
        if (prec .le. zero) then
            call utmess('F', 'MODELISA3_98')
        end if
    end if
!
! --- RECUPERATION DU NOMBRE DE NOEUDS DU MAILLAGE :
!     ---------------------------------------------
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoe)
!
! --- ALLOCATION DU VECTEUR DES NOMS DES NOEUDS  APPARTENANT
! --- AU PLAN OU A LA DROITE :
!     ----------------------
    call wkvect(lisnoe, 'V V I', nbnoe, idlino)
!
    call cgnop0(nbnoe, vale, x0, vecnor, prec, &
                nbno, zi(idlino))
!
    call jedema()
!.============================ FIN DE LA ROUTINE ======================
end subroutine
