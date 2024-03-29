! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine te0022(option, nomte)
    implicit none
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DES CONTRAINTES AUX POINTS DE GAUSS
!          ELEMENTS ISOPARAMETRIQUES 2D et 3D
!
!          OPTION : 'SIEF_ELGA'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/ortrep.h"
#include "asterfort/sigvmc.h"
!
    integer :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano
    integer :: idim
    integer :: i, icont, idepl, igeom, imate, nbsig
!
    real(kind=8) :: sigma(162), repere(7), instan, nharm
    real(kind=8) :: bary(3)
    real(kind=8) :: zero
!
!-----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! - NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!   -----------------------------------------
    nbsig = nbsigm()
!
! - INITIALISATIONS :
!   -----------------
    zero = 0.0d0
    instan = r8vide()
    nharm = zero
!
    do i = 1, nbsig*npg
        sigma(i) = zero
    end do
!
! - RECUPERATION DES COORDONNEES DES CONNECTIVITES
!   ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! - RECUPERATION DU MATERIAU
!   ------------------------
    call jevech('PMATERC', 'L', imate)
!
! - RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE
!   ------------------------------------------------------------
!   COORDONNEES DU BARYCENTRE ( POUR LE REPRE CYLINDRIQUE )
!
    bary(1) = 0.d0
    bary(2) = 0.d0
    bary(3) = 0.d0
    do i = 1, nno
        do idim = 1, ndim
            bary(idim) = bary(idim)+zr(igeom+idim+ndim*(i-1)-1)/nno
        end do
    end do
    call ortrep(ndim, bary, repere)
!
! ---- RECUPERATION DU CHAMP DE DEPLACEMENT SUR L'ELEMENT
!      --------------------------------------------------
    call jevech('PDEPLAR', 'L', idepl)
!
    call sigvmc('RIGI', nno, ndim, nbsig, npg, &
                ipoids, ivf, idfde, zr(igeom), zr(idepl), &
                instan, repere, zi(imate), nharm, sigma)
!
!
! ---- RECUPERATION ET AFFECTATION DU VECTEUR EN SORTIE
! ---- AVEC LE VECTEUR DES CONTRAINTES AUX POINTS D'INTEGRATION
!      --------------------------------------------------------
    call jevech('PCONTRR', 'E', icont)
!
    do i = 1, nbsig*npg
        zr(icont+i-1) = sigma(i)
    end do
!
end subroutine
