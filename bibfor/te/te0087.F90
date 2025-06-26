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
subroutine te0087(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epsvmc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!
!     BUT: CALCUL DES DEFORMATIONS AUX POINTS D'INTEGRATION
!          DES ELEMENTS ISOPARAMETRIQUES 2D
!
!          OPTIONS : 'EPSI_ELGA'
!                    'EPSG_ELGA'
!                    'EPME_ELGA'
!                    'EPMG_ELGA'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: nbsig, nbsig1, nbsig2, ndim, nno, i
    integer(kind=8) :: nnos, npg, ipoids, ivf, idfde
    integer(kind=8) :: igau, isig, igeom, idepl, iret
    integer(kind=8) :: itemps, idefo
    integer(kind=8) :: jgano
!
    real(kind=8) :: epsm(54), angl_naut(3)
    real(kind=8) :: nharm, instan, zero
!
    character(len=4) :: fami
!
!
! DEB ------------------------------------------------------------------
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig1 = nbsigm()
    nbsig2 = 6
    nbsig = nbsig1
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    instan = r8vide()
    nharm = zero
!
    do i = 1, nbsig2*npg
        epsm(i) = zero
    end do
!
!
! ---- RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!      ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! ---- RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE :
!      ------------------------------------------------------------
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! ---- RECUPERATION DU CHAMP DE DEPLACEMENT SUR L'ELEMENT :
!      --------------------------------------------------
    call jevech('PDEPLAR', 'L', idepl)
!
! ---- RECUPERATION DE L'INSTANT DE CALCUL :
!      -----------------------------------
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) then
        instan = zr(itemps)
    end if
!
! ---- RECUPERATION DU VECTEUR DES DEFORMATIONS EN SORTIE :
!      --------------------------------------------------
    call jevech('PDEFOPG', 'E', idefo)
!
! ---- CALCUL DES DEFORMATIONS MECANIQUES AUX POINTS D'INTEGRATION
! ---- DE L'ELEMENT , I.E. SI ON NOTE EPSI_MECA = B*U
! ---- ON CALCULE SIMPLEMENT EPSI_MECA POUR LES OPTIONS EPSI ET EPSG
! ----                    ET EPSI_MECA - EPSI_THERMIQUES POUR LES
! ----                    OPTIONS EPME ET EPMG :
!      ---------------------------------------
    call epsvmc(fami, nno, ndim, nbsig1, npg, &
                ipoids, ivf, idfde, zr(igeom), zr(idepl), &
                instan, angl_naut, nharm, option, epsm)
!
!      --------------------
! ---- AFFECTATION DU VECTEUR EN SORTIE AVEC LES DEFORMATIONS AUX
! ---- POINTS D'INTEGRATION :
!      --------------------
    do igau = 1, npg
        do isig = 1, nbsig
            zr(idefo+nbsig*(igau-1)+isig-1) = epsm(nbsig*(igau-1)+isig)
        end do
    end do
!
end subroutine
