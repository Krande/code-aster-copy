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
subroutine eps2mc(nno, ndim, nbsig, npg, ipoids, &
                  ivf, idfde, xyz, depl, eps2)
!.======================================================================
    implicit none
!
!      EPS2MC   -- CALCUL DES  DEFORMATIONS DU SECOND ORDRE AUX
!                  POINTS D'INTEGRATION POUR LES ELEMENTS
!                  ISOPARAMETRIQUES
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NNO            IN     I        NOMBRE DE NOEUDS DE L'ELEMENT
!    NDIM           IN     I        DIMENSION DE L'ELEMENT (2 OU 3)
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE
!                                   A L'ELEMENT
!    NPG            IN     I        NOMBRE DE POINTS D'INTEGRATION
!                                   DE L'ELEMENT
!    IVF            IN     I        POINTEUR FONCTIONS DE FORME
!    IPOIDS         IN     I        POINTEUR POIDS D'INTEGRATION
!    IDFDE          IN     I        PT DERIVEES DES FONCTIONS DE FORME
!    XYZ(1)         IN     R        COORDONNEES DES CONNECTIVITES
!    DEPL(1)        IN     R        VECTEUR DES DEPLACEMENTS SUR
!                                   L'ELEMENT
!    EPS2(1)        OUT    R        DEFORMATIONS DU SECOND ORDRE
!                                   AUX POINTS D'INTEGRATION
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
    real(kind=8) :: xyz(1), depl(1), eps2(1)
! -----  VARIABLES LOCALES
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27)
    real(kind=8) :: jacob
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     -----------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idecno, idfde, igau, ipoids, ivf, k
    integer(kind=8) :: nbsig, ndim, nno, npg
    real(kind=8) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx
    real(kind=8) :: dwdy, dwdz, dx, rayon, undemi, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    undemi = 0.5d0
!
    do i = 1, nbsig*npg
        eps2(i) = zero
    end do
!
! --- CALCUL DES DEFORMATIONS DU SECOND ORDRE AUX POINTS D'INTEGRATION
! ---  BOUCLE SUR LES POINTS D'INTEGRATION :
!      -----------------------------------
    do igau = 1, npg
!
        dx = zero
        rayon = zero
        dudx = zero
        dudy = zero
        dudz = zero
        dvdx = zero
        dvdy = zero
        dvdz = zero
        dwdx = zero
        dwdy = zero
        dwdz = zero
!
!       -------------
! ----  CAS MASSIF 3D
!       -------------
        if (lteatt('DIM_TOPO_MAILLE', '3')) then
!
! ----    CALCUL DES DERIVEES DES FONCTIONS DE FORME SUR L'ELEMENT
! ----    REEL ET DU PRODUIT JACOBIEN*POIDS (DANS JACOB) :
!         ----------------------------------------------
            call dfdm3d(nno, igau, ipoids, idfde, xyz, &
                        jacob, dfdx, dfdy, dfdz)
!
! ----    CALCUL DES DERIVEES DES DEPLACEMENTS :
!         ------------------------------------
            do i = 1, nno
!
                dudx = dudx+dfdx(i)*depl((i-1)*ndim+1)
                dudy = dudy+dfdy(i)*depl((i-1)*ndim+1)
                dudz = dudz+dfdz(i)*depl((i-1)*ndim+1)
!
                dvdx = dvdx+dfdx(i)*depl((i-1)*ndim+2)
                dvdy = dvdy+dfdy(i)*depl((i-1)*ndim+2)
                dvdz = dvdz+dfdz(i)*depl((i-1)*ndim+2)
!
                dwdx = dwdx+dfdx(i)*depl((i-1)*ndim+3)
                dwdy = dwdy+dfdy(i)*depl((i-1)*ndim+3)
                dwdz = dwdz+dfdz(i)*depl((i-1)*ndim+3)
!
            end do
!
! ----    DEFORMATIONS DU SECOND ORDRE :
!         ----------------------------
            eps2(nbsig*(igau-1)+1) = undemi*(dudx*dudx+dvdx*dvdx+dwdx*dwdx)
            eps2(nbsig*(igau-1)+2) = undemi*(dudy*dudy+dvdy*dvdy+dwdy*dwdy)
            eps2(nbsig*(igau-1)+3) = undemi*(dudz*dudz+dvdz*dvdz+dwdz*dwdz)
!
            eps2(nbsig*(igau-1)+4) = undemi*(dudx*dudy+dvdx*dvdy+dwdx*dwdy)
            eps2(nbsig*(igau-1)+5) = undemi*(dudx*dudz+dvdx*dvdz+dwdx*dwdz)
            eps2(nbsig*(igau-1)+6) = undemi*(dudy*dudz+dvdy*dvdz+dwdy*dwdz)
!
!       ------------------------------------------------------------
! ----  CAS MASSIF 2D CONTRAINTES PLANES, DEFORMATIONS PLANES ET AXI
!       ------------------------------------------------------------
        elseif (lteatt('C_PLAN', 'OUI') .or. lteatt('D_PLAN', &
                                                    'OUI') .or. lteatt('AXIS', 'OUI')) then
!
            k = (igau-1)*nno
!
! ----    CALCUL DES DERIVEES DES FONCTIONS DE FORME SUR L'ELEMENT
! ----    REEL ET DU PRODUIT JACOBIEN*POIDS (DANS JACOB) :
!         ----------------------------------------------
            call dfdm2d(nno, igau, ipoids, idfde, xyz, &
                        jacob, dfdx, dfdy)
!
! ----    CALCUL DES DERIVEES DES DEPLACEMENTS :
!         ------------------------------------
            do i = 1, nno
!
                dudx = dudx+dfdx(i)*depl((i-1)*ndim+1)
                dudy = dudy+dfdy(i)*depl((i-1)*ndim+1)
!
                dvdx = dvdx+dfdx(i)*depl((i-1)*ndim+2)
                dvdy = dvdy+dfdy(i)*depl((i-1)*ndim+2)
!
                if (lteatt('AXIS', 'OUI')) then
                    idecno = 2*(i-1)
                    rayon = rayon+zr(ivf+i+k-1)*xyz(1+idecno)
                    dx = dx+zr(ivf+i+k-1)*depl(1+idecno)
                end if
            end do
!
! ----    DEFORMATIONS DU SECOND ORDRE :
!         ----------------------------
            eps2(nbsig*(igau-1)+1) = undemi*(dudx*dudx+dvdx*dvdx)
            eps2(nbsig*(igau-1)+2) = undemi*(dudy*dudy+dvdy*dvdy)
            eps2(nbsig*(igau-1)+3) = zero
!
            if (lteatt('AXIS', 'OUI')) then
                eps2(nbsig*(igau-1)+3) = undemi*dx*dx/rayon/rayon
            end if
!
            eps2(nbsig*(igau-1)+4) = undemi*(dudx*dudy+dvdx*dvdy)
        else
            call utmess('F', 'ELEMENTS_11')
        end if
!
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
