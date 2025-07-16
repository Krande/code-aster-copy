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
subroutine te0198(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/sigtmc.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_MECA_TEMP_R  '
!                          ELEMENTS FOURIER
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    character(len=4) :: fami
    real(kind=8) :: bsigma(81), sigth(162), angl_naut(3)
    real(kind=8) :: nharm, instan
    integer(kind=8) :: ndim, nno, nnos, npg1, ipoids, ivf, idfde, jgano
    integer(kind=8) :: dimmod
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, igeom, iharmo, imate, iret, itemps, ivectu
    integer(kind=8) :: nbsig, nh
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    instan = r8vide()
    nharm = zero
    dimmod = 3
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
    do i = 1, nbsig*npg1
        sigth(i) = zero
    end do
!
    do i = 1, dimmod*nno
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
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! ---- RECUPERATION DE L'INSTANT
!      -------------------------
    call tecach('NNO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) instan = zr(itemps)
!
! ---- RECUPERATION  DU NUMERO D'HARMONIQUE
!      ------------------------------------
    call jevech('PHARMON', 'L', iharmo)
    nh = zi(iharmo)
    nharm = dble(nh)
!
! ---- CALCUL DES CONTRAINTES THERMIQUES AUX POINTS D'INTEGRATION
! ---- DE L'ELEMENT :
!      ------------
    call sigtmc(fami, dimmod, nbsig, npg1, &
                instan, zi(imate), angl_naut, &
                option, sigth)
!
! ---- CALCUL DU VECTEUR DES FORCES D'ORIGINE THERMIQUE/HYDRIQUE
! ---- OU DE SECHAGE (BT*SIGTH)
!      ----------------------------------------------------------
    call bsigmc(nno, dimmod, nbsig, npg1, ipoids, &
                ivf, idfde, zr(igeom), nharm, sigth, &
                bsigma)
!
! ---- RECUPERATION ET AFFECTATION DU VECTEUR EN SORTIE AVEC LE
! ---- VECTEUR DES FORCES D'ORIGINE THERMIQUE
!      -------------------------------------
    call jevech('PVECTUR', 'E', ivectu)
!
    do i = 1, dimmod*nno
        zr(ivectu+i-1) = bsigma(i)
    end do
!
! FIN ------------------------------------------------------------------
end subroutine
