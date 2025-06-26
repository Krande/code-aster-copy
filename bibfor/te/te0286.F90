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
subroutine te0286(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ethdst.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/simtep.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
!.......................................................................
! FONCTION REALISEE:
!
!      CALCUL DE L'ENERGIE POTENTIELLE THERMOELASTIQUE A L'EQUILIBRE
!      ELEMENTS ISOPARAMETRIQUES 2D
!
!      OPTION : 'EPOT_ELEM'
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    character(len=4) :: fami
    real(kind=8) :: sigma(162), bsigma(81), angl_naut(3)
    real(kind=8) :: instan, nharm
!
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idepl, idfde, iener, igeom, iharmo, imate
    integer(kind=8) :: ipoids, iret, ivf, jgano, nbsig, ndim, ndim2
    integer(kind=8) :: nh, nno, nnos, npg
    real(kind=8) :: enthth, epot, undemi, zero
!-----------------------------------------------------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! --- INITIALISATIONS :
!     -----------------
    nh = 0
    zero = 0.0d0
    undemi = 0.5d0
    instan = r8vide()
    nharm = zero
    ndim2 = 2
    if (lteatt('FOURIER', 'OUI')) then
        ndim = 3
    end if
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
    do i = 1, nbsig*npg
        sigma(i) = zero
    end do
!
    do i = 1, ndim*nno
        bsigma(i) = zero
    end do
!
! ---- RECUPERATION DES COORDONNEES DES CONNECTIVITES
!      ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! ---- RECUPERATION DU MATERIAU
!      ------------------------
    call jevech('PMATERC', 'L', imate)
!
! ---- RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE
!      ------------------------------------------------------------
    call getElemOrientation(ndim2, nno, igeom, angl_naut)
!
! ---- RECUPERATION DU CHAMP DE DEPLACEMENT SUR L'ELEMENT
!      --------------------------------------------------
    call jevech('PDEPLAR', 'L', idepl)
!
! ---- RECUPERATION  DU NUMERO D'HARMONIQUE
!      ------------------------------------
    call tecach('NNO', 'PHARMON', 'L', iret, iad=iharmo)
    if (iharmo .ne. 0) then
        nh = zi(iharmo)
        nharm = dble(nh)
    end if
!
! ---- CALCUL DES CONTRAINTES 'VRAIES' SUR L'ELEMENT
! ---- (I.E.  1/2*SIGMA_MECA - SIGMA_THERMIQUES)
!      ------------------------------------
    call simtep(fami, nno, ndim, nbsig, npg, &
                ipoids, ivf, idfde, zr(igeom), zr(idepl), &
                instan, angl_naut, zi(imate), nharm, sigma)
!
! ---- CALCUL DU VECTEUR DES FORCES INTERNES (BT*SIGMA)
!      -----------------------------------------------
    call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                ivf, idfde, zr(igeom), nharm, sigma, &
                bsigma)
!
! ---- CALCUL DU TERME EPSTH_T*D*EPSTH
!      -------------------------------
    call ethdst(fami, nno, ndim, nbsig, npg, &
                ipoids, ivf, idfde, zr(igeom), zr(idepl), &
                instan, angl_naut, zi(imate), option, enthth)
!
! ---- CALCUL DE L'ENERGIE POTENTIELLE :
! ----        1/2*UT*K*U - UT*FTH + 1/2*EPSTHT*D*EPSTH :
!             ----------------------------------------
!
    epot = zero
!
    do i = 1, ndim*nno
        epot = epot+bsigma(i)*zr(idepl+i-1)
    end do
!
    epot = epot+undemi*enthth
!
! ---- RECUPERATION ET AFFECTATION DU REEL EN SORTIE
! ---- AVEC L'ENERGIE DE DEFORMATION
!      -----------------------------
    call jevech('PENERDR', 'E', iener)
!
    zr(iener) = epot
!
end subroutine
