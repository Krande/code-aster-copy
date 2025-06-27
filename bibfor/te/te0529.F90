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
subroutine te0529(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epstmc.h"
#include "asterfort/jevech.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
!
!
!     BUT: CALCUL DES DEFORMATIONS LIEES AUX VARIABLES DE COMMANDE
!          AUX POINTS D'INTEGRATION DES ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTIONS : 'EPVC_ELGA'
!    CINQ COMPOSANTES :
!    EPTHER_L = DILATATION THERMIQUE (LONGI)   : ALPHA_L*(T-TREF)
!    EPTHER_T = DILATATION THERMIQUE (TRANSV)   : ALPHA_T*(T-TREF)
!    EPTHER_N = DILATATION THERMIQUE (NORMLALE)   : ALPHA_N*(T-TREF)
!    EPSECH = RETRAIT DE DESSICCATION : -K_DESSIC(SREF-SECH)
!    EPHYDR = RETRAIT ENDOGENE        : -B_ENDOGE*HYDR
!    EPPTOT = RETRAIT DU A LA PRESSION DE FLUIDE EN THM CHAINEE :
!             -BIOT*PTOT
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    integer(kind=8) :: jgano, ndim, nno, i, nnos, npg, ipoids, ivf, idfde, igau, isig
    integer(kind=8) :: igeom, itemps, idefo, imate, iret, nbcmp
    real(kind=8) :: epvc(162), angl_naut(3)
    real(kind=8) :: instan, epsse(6), epsth(6), epshy(6), epspt(6)
    character(len=4) :: fami
    character(len=16) :: optio2
! DEB ------------------------------------------------------------------
!
!    NOMBRE DE COMPOSANTES  A  CALCULER
    nbcmp = 6
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! ---- RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!      ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! ---- RECUPERATION DU MATERIAU :
!      ----------------------------------------------
    call tecach('NNO', 'PMATERC', 'L', iret, iad=imate)
!
! --- RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE :
!     ------------------------------------------------------------
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! ---- RECUPERATION DE L'INSTANT DE CALCUL :
!      -----------------------------------
    call tecach('NNO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) then
        instan = zr(itemps)
    else
        instan = r8vide()
    end if
!
!     -----------------
! ---- RECUPERATION DU VECTEUR DES DEFORMATIONS EN SORTIE :
!      --------------------------------------------------
    call jevech('PDEFOPG', 'E', idefo)
    call r8inir(135, 0.d0, epvc, 1)
!
!
    do igau = 1, npg
!
!      CALCUL AU POINT DE GAUSS DE LA TEMPERATURE
! ------------------------------------------
!
        optio2 = 'EPVC_ELGA_TEMP'
!
        call epstmc(fami, ndim, instan, '+', igau, &
                    1, angl_naut, zi(imate), optio2, &
                    epsth)
        optio2 = 'EPVC_ELGA_SECH'
        call epstmc(fami, ndim, instan, '+', igau, &
                    1, angl_naut, zi(imate), optio2, &
                    epsse)
        optio2 = 'EPVC_ELGA_HYDR'
        call epstmc(fami, ndim, instan, '+', igau, &
                    1, angl_naut, zi(imate), optio2, &
                    epshy)
        optio2 = 'EPVC_ELGA_PTOT'
        call epstmc(fami, ndim, instan, '+', igau, &
                    1, angl_naut, zi(imate), optio2, &
                    epspt)
        do i = 1, 3
            epvc(i+nbcmp*(igau-1)) = epsth(i)
        end do
        epvc(4+nbcmp*(igau-1)) = epsse(1)
        epvc(5+nbcmp*(igau-1)) = epshy(1)
        epvc(6+nbcmp*(igau-1)) = epspt(1)
!
    end do
!
!         --------------------
! ---- AFFECTATION DU VECTEUR EN SORTIE AVEC LES DEFORMATIONS AUX
! ---- POINTS D'INTEGRATION :
!      --------------------
    do igau = 1, npg
        do isig = 1, nbcmp
            zr(idefo+nbcmp*(igau-1)+isig-1) = epvc(nbcmp*(igau-1)+isig)
        end do
    end do
!
!
end subroutine
