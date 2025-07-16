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
subroutine te0426(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/rcvarc.h"
#include "asterfort/sigimc.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_MECA_EPSA_R  '
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    character(len=4) :: fami
!
    real(kind=8) :: sigi(162), epsi(162), bsigma(81), angl_naut(3)
    real(kind=8) :: instan, nharm
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idfde, igau, igeom, imate, ipoids, iret
    integer(kind=8) :: itemps, ivectu, ivf, jgano, nbsig, ndim, nno
    integer(kind=8) :: nnos, npg1
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    instan = r8vide()
    nharm = zero
!
    do i = 1, nbsig*npg1
        epsi(i) = zero
        sigi(i) = zero
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
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! ---- RECUPERATION DE L'INSTANT
!      -------------------------
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) instan = zr(itemps)
!
!
! ---- CONSTRUCTION DU VECTEUR DES DEFORMATIONS ANELASTIQUES DEFINIES
! ---- AUX POINTS D'INTEGRATION A PARTIR DES DONNEES UTILISATEUR
! ---- + MISE AU FORMAT DES TERMES EXTRA-DIAGONAUX (COHERENT AVEC DMATMC)
!      --------------------------------------------------------------
!
    do igau = 1, npg1
        call rcvarc(' ', 'EPSAXX', '+', 'RIGI', igau, &
                    1, epsi(nbsig*(igau-1)+1), iret)
        if (iret .eq. 1) epsi(nbsig*(igau-1)+1) = 0.d0
!
        call rcvarc(' ', 'EPSAYY', '+', 'RIGI', igau, &
                    1, epsi(nbsig*(igau-1)+2), iret)
        if (iret .eq. 1) epsi(nbsig*(igau-1)+2) = 0.d0
!
        call rcvarc(' ', 'EPSAZZ', '+', 'RIGI', igau, &
                    1, epsi(nbsig*(igau-1)+3), iret)
        if (iret .eq. 1) epsi(nbsig*(igau-1)+3) = 0.d0
!
        call rcvarc(' ', 'EPSAXY', '+', 'RIGI', igau, &
                    1, epsi(nbsig*(igau-1)+4), iret)
        if (iret .eq. 1) epsi(nbsig*(igau-1)+4) = 0.d0
        epsi(nbsig*(igau-1)+4) = 2.0*epsi(nbsig*(igau-1)+4)
!
        call rcvarc(' ', 'EPSAXZ', '+', 'RIGI', igau, &
                    1, epsi(nbsig*(igau-1)+5), iret)
        if (iret .eq. 1) epsi(nbsig*(igau-1)+5) = 0.d0
        epsi(nbsig*(igau-1)+5) = 2.0*epsi(nbsig*(igau-1)+5)
!
        call rcvarc(' ', 'EPSAYZ', '+', 'RIGI', igau, &
                    1, epsi(nbsig*(igau-1)+6), iret)
        if (iret .eq. 1) epsi(nbsig*(igau-1)+6) = 0.d0
        epsi(nbsig*(igau-1)+6) = 2.0*epsi(nbsig*(igau-1)+6)
    end do
!
! ---- CALCUL DU VECTEUR DES CONTRAINTES ANELASTIQUES AUX POINTS
! ---- D'INTEGRATION
!      -------------
    call sigimc(fami, nno, ndim, nbsig, npg1, &
                instan, zi(imate), angl_naut, &
                epsi, sigi)
!
! ---- CALCUL DU VECTEUR DES FORCES DUES AUX CONTRAINTES ANELASTIQUES
! ---- (I.E. BT*SIG_ANELASTIQUES)
!      ----------------------
    call bsigmc(nno, ndim, nbsig, npg1, ipoids, &
                ivf, idfde, zr(igeom), nharm, sigi, &
                bsigma)
!
! ---- RECUPERATION ET AFFECTATION DU VECTEUR EN SORTIE AVEC LE
! ---- VECTEUR DES FORCES DUES AUX CONTRAINTES ANELASTIQUES
!      -------------------------------------------------
    call jevech('PVECTUR', 'E', ivectu)
!
    do i = 1, ndim*nno
        zr(ivectu+i-1) = bsigma(i)
    end do
!
! FIN ------------------------------------------------------------------
end subroutine
