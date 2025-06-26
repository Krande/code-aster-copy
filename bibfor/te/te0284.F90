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

subroutine te0284(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epsimc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/sigimc.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES EN 2D
!                      OPTION : 'CHAR_MECA_EPSI_R  ','CHAR_MECA_EPSI_F '
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    character(len=4) :: fami
    real(kind=8) :: sigi(162), epsi(162), bsigma(81), angl_naut(3)
    real(kind=8) :: instan, nharm, xyz(81)
    integer(kind=8) :: dimcoo, idim
!
!
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idfde, igeom, iharmo, imate, ipoids, iret
    integer(kind=8) :: itemps, ivectu, ivf, jgano, nbsig, ndim, nno
    integer(kind=8) :: nnos, npg
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    dimcoo = ndim
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    instan = r8vide()
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
    if (nbsig .eq. 6) ndim = 3
!
    do i = 1, nbsig*npg
        epsi(i) = zero
        sigi(i) = zero
    end do
!
    do i = 1, ndim*nno
        bsigma(i) = zero
    end do
!
! ---- RECUPERATION DE L'HARMONIQUE DE FOURIER
!      ---------------------------------------
    call tecach('NNO', 'PHARMON', 'L', iret, iad=iharmo)
    if (iharmo .eq. 0) then
        nharm = zero
    else
        nharm = dble(zi(iharmo))
    end if
!
! ---- RECUPERATION DES COORDONNEES DES CONNECTIVITES
!      ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
    if (ndim .eq. dimcoo) then
        do i = 1, ndim*nno
            xyz(i) = zr(igeom+i-1)
        end do
    else
        do i = 1, nno
            do idim = 1, ndim
                if (idim .le. dimcoo) then
                    xyz(idim+ndim*(i-1)) = zr(igeom-1+idim+dimcoo*(i-1))
                else
                    xyz(idim+ndim*(i-1)) = 0.d0
                end if
            end do
        end do
    end if
!
! ---- RECUPERATION DU MATERIAU
!      ------------------------
    call jevech('PMATERC', 'L', imate)
!
! ---- RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE
!      ------------------------------------------------------------
    call getElemOrientation(dimcoo, nno, igeom, angl_naut)
!
!
! ---- RECUPERATION DE L'INSTANT
!      -------------------------
    call tecach('NNO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) instan = zr(itemps)
!
! ---- CONSTRUCTION DU VECTEUR DES DEFORMATIONS INITIALES DEFINIES AUX
! ---- POINTS D'INTEGRATION A PARTIR DES DONNEES UTILISATEUR
!      -----------------------------------------------------
    call epsimc(option, zr(igeom), nno, npg, ndim, &
                nbsig, zr(ivf), epsi)
!
! ---- CALCUL DU VECTEUR DES CONTRAINTES INITIALES AUX POINTS
! ---- D'INTEGRATION
!      -------------
    call sigimc(fami, nno, ndim, nbsig, npg, &
                instan, zi(imate), angl_naut, &
                epsi, sigi)
!
! ---- CALCUL DU VECTEUR DES FORCES DUES AUX CONTRAINTES INITIALES
! ---- (I.E. BT*SIG_INITIALES)
!      ----------------------
    call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                ivf, idfde, zr(igeom), nharm, sigi, &
                bsigma)
!
! ---- RECUPERATION ET AFFECTATION DU VECTEUR EN SORTIE AVEC LE
! ---- VECTEUR DES FORCES DUES AUX CONTRAINTES INITIALES
!      -------------------------------------------------
    call jevech('PVECTUR', 'E', ivectu)
!
    do i = 1, ndim*nno
        zr(ivectu+i-1) = bsigma(i)
    end do
!
! FIN ------------------------------------------------------------------
end subroutine
