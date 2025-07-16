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
subroutine ethdst(fami, nno, ndim, nbsig, npg, &
                  ipoids, ivf, idfde, xyz, depl, &
                  instan, angl_naut, mater, option, enthth)
    implicit none
!
!      ETHDST   -- CALCUL DU TERME EPSTHT*D*EPSTH RENTRANT
!                  DANS LE CALCUL DE L'ENERGIE POTENTIELLE
!                  (I.E.  1/2*UT*K*U - UT*FTH + 1/2*EPSTHT*D*EPSTH)
!                  POUR LES ELEMENTS ISOPARAMETRIQUES
!
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NNO            IN     I        NOMBRE DE NOEUDS DE L'ELEMENT
!    NDIM           IN     I        DIMENSION DE L'ELEMENT (2 OU 3)
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE
!                                   A L'ELEMENT
!    NPG            IN     I        NOMBRE DE POINTS D'INTEGRATION
!                                   DE L'ELEMENT
!    IPOIDS         IN     I        POIDS D'INTEGRATION
!    IVF            IN     I        FONCTIONS DE FORME
!    IDFDE          IN     I        DERIVEES DES FONCTIONS DE FORME
!    XYZ(1)         IN     R        COORDONNEES DES CONNECTIVITES
!    DEPL(1)        IN     R        VECTEUR DES DEPLACEMENTS SUR
!                                   L'ELEMENT
!    INSTAN         IN     R        INSTANT DE CALCUL
!    ANGL_NAUT(3)   IN     R        ANGLES NAUTIQUES DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    MATER          IN     I        MATERIAU
!    OPTION         IN     K16      OPTION DE CALCUL
!    ENTHTH         OUT    R        SOMME(EPSTH_T*D*EPSTH)
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/epthmc.h"
#include "asterfort/lteatt.h"
#include "asterfort/sigtmc.h"
!
    integer(kind=8) :: ipoids, ivf, idfde
    character(len=16) :: option
    character(len=*) :: fami
    real(kind=8) :: xyz(*), depl(*), angl_naut(3)
    real(kind=8) :: instan, enthth
! -----  VARIABLES LOCALES
    integer(kind=8) :: i, mater, nbsig, ndim, nno, npg, k, igau
    character(len=16) :: k16bid
    real(kind=8) :: sigth(162), zero
    real(kind=8) :: rayon
    real(kind=8) :: epsith(162), enthpg, dfdx(27), dfdy(27), dfdz(27)
    real(kind=8) :: poidi
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    k16bid = ' '
    enthth = zero
!
! --- CALCUL DES CONTRAINTES MECANIQUES AUX POINTS D'INTEGRATION
!      ---------------------------------------------------------
    call epthmc(fami, nno, ndim, nbsig, npg, &
                zr(ivf), angl_naut, instan, mater, &
                option, epsith)
!
! --- CALCUL DES CONTRAINTES THERMIQUES AUX POINTS D'INTEGRATION
!      ---------------------------------------------------------
    call sigtmc(fami, ndim, nbsig, npg, &
                instan, mater, angl_naut, &
                k16bid, sigth)
!
! --- CALCUL DES CONTRAINTES TOTALES AUX POINTS D'INTEGRATION
!      ---------------------------------------------------------
    do igau = 1, npg
        enthpg = 0.d0
! ----  CALCUL DU JACOBIEN*POIDS - CAS MASSIF 3D
!
        if (lteatt('DIM_TOPO_MAILLE', '3')) then
            call dfdm3d(nno, igau, ipoids, idfde, xyz, &
                        poidi, dfdx, dfdy, dfdz)
! ----  CALCUL DU JACOBIEN*POIDS - CAS MASSIF 2D
        else
            k = (igau-1)*nno
            call dfdm2d(nno, igau, ipoids, idfde, xyz, &
                        poidi, dfdx, dfdy)
            if (lteatt('AXIS', 'OUI')) then
                rayon = 0.d0
                do i = 1, nno
                    rayon = rayon+zr(ivf+k-1+i)*xyz(2*(i-1)+1)
                end do
                poidi = poidi*rayon
            end if
        end if
        do i = 1, nbsig
            enthpg = enthpg+epsith(i+nbsig*(igau-1))*sigth(i+nbsig*(igau-1))
        end do
        enthth = enthth+(enthpg*poidi)
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
