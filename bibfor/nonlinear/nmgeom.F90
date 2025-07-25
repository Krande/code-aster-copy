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
subroutine nmgeom(ndim, nno, axi, grand, geom, &
                  kpg, ipoids, ivf, idfde, depl, &
                  ldfdi, poids, dfdi, f, eps, &
                  r)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
    aster_logical :: axi, grand
    integer(kind=8) :: ndim, nno, kpg
    real(kind=8) :: geom(ndim, nno), dfdi(nno, ndim), depl(ndim, nno)
    real(kind=8) :: poids, f(3, 3), eps(6), r
    aster_logical :: ldfdi
!
!.......................................................................
!
!     BUT:  CALCUL DES ELEMENTS CINEMATIQUES (MATRICES F ET E, RAYON R)
!           EN UN POINT DE GAUSS (EVENTUELLEMENT EN GRANDES TRANSFORM.)
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  AXI     : INDICATEUR SI AXISYMETRIQUE
! IN  GRAND   : INDICATEUR SI GRANDES TRANSFORMATIONS
! IN  GEOM    : COORDONEES DES NOEUDS
! IN  KPG     : NUMERO DU POINT DE GAUSS (POUR L'ACCES AUX FCT. FORMES)
! IN  IPOIDS  : POIDS DU POINT DE GAUSS DE L'ELEMENT DE REFERENCE
! IN  IVF     : VALEUR DES FONCTIONS DE FORME (EN AXISYMETRIQUE)
! IN  IDFDE   : DERIVEE DES FONCTIONS DE FORME DE REFERENCE
! IN  DEPL    : DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  DEPL    : DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  LDFDI   : VEUT-ON CALCULER DFDI ET POIDS
! OUT POIDS   : "POIDS" DU POINT DE GAUSS
! OUT DFDI    : DERIVEE DES FONCTIONS DE FORME
! OUT F       : GRADIENT DE LA TRANSFORMATION
! OUT EPS     : DEFORMATIONS
! OUT R       : DISTANCE DU POINT DE GAUSS A L'AXE (EN AXISYMETRIQUE)
!......................................................................
! REMARQUE CONCERNANT L'ARGUMENT LDFDI :
!  NMGEOM EST PARFOIS APPELE 2 FOIS DE SUITE AVEC U ET DELTA_U (PAR
!  EXEMPLE DANS NMPL3D). COMME LE CALCUL DE DFDM3D EST COUTEUX ET QU'IL
!  EST INDEPENDANT DE U, ON PEUT ECONOMISER LE 2EME CALCUL EN UTILISANT
!  L'ARGUMENT LDFDI : 1ER APPEL .TRUE. ; 2EME APPEL .FALSE.
!
!
    aster_logical :: tridim
    integer(kind=8) :: i, j, k, n
    real(kind=8) :: grad(3, 3), epstab(3, 3), ur, tmp
    real(kind=8) :: rac2, kron(3, 3)
!-----------------------------------------------------------------------
    integer(kind=8) :: idfde, ipoids, ivf
!-----------------------------------------------------------------------
    data kron/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
    rac2 = sqrt(2.d0)
    tridim = (ndim .eq. 3)
!
! - CALCUL DES DERIVEES DES FONCTIONS DE FORME ET JACOBIEN
    if (ldfdi) then
        if (tridim) then
            call dfdm3d(nno, kpg, ipoids, idfde, geom, &
                        poids, dfdi(1, 1), dfdi(1, 2), dfdi(1, 3))
        else
            call dfdm2d(nno, kpg, ipoids, idfde, geom, &
                        poids, dfdi(1, 1), dfdi(1, 2))
        end if
    end if
!
!
! - CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE) ET DU DEPL. RADIAL
    if (axi) then
        r = 0.d0
        ur = 0.d0
        do n = 1, nno
            r = r+zr(ivf-1+n+(kpg-1)*nno)*geom(1, n)
            ur = ur+zr(ivf-1+n+(kpg-1)*nno)*depl(1, n)
        end do
        if (ldfdi) poids = poids*r
    end if
!
! - CALCUL DES GRADIENT : GRAD(U) ET F
!
    do i = 1, 3
        do j = 1, 3
            f(i, j) = kron(i, j)
            grad(i, j) = 0.d0
        end do
    end do
!
    if (tridim) then
        do n = 1, nno
            do i = 1, 3
                do j = 1, 3
                    grad(i, j) = grad(i, j)+dfdi(n, j)*depl(i, n)
                end do
            end do
        end do
    else
        do n = 1, nno
            do i = 1, 2
                do j = 1, 2
                    grad(i, j) = grad(i, j)+dfdi(n, j)*depl(i, n)
                end do
            end do
        end do
    end if
!
    if (grand) then
        do i = 1, 3
            do j = 1, 3
                f(i, j) = f(i, j)+grad(i, j)
            end do
        end do
        if (axi) f(3, 3) = 1.d0+ur/r
    end if
!
! - CALCUL DES DEFORMATIONS : E
!
    do i = 1, ndim
        do j = 1, i
            tmp = grad(i, j)+grad(j, i)
!
            if (grand) then
                do k = 1, ndim
                    tmp = tmp+grad(k, i)*grad(k, j)
                end do
            end if
!
            epstab(i, j) = 0.5d0*tmp
!
        end do
    end do
!
    eps(1) = epstab(1, 1)
    eps(2) = epstab(2, 2)
    eps(3) = 0.d0
    eps(4) = epstab(2, 1)*rac2
    eps(5) = 0.d0
    eps(6) = 0.d0
!
    if (tridim) then
        eps(3) = epstab(3, 3)
        eps(5) = epstab(3, 1)*rac2
        eps(6) = epstab(3, 2)*rac2
    else if (axi) then
        eps(3) = ur/r
        if (grand) eps(3) = eps(3)+0.5d0*ur*ur/(r*r)
    else
        eps(3) = 0.d0
    end if
!
end subroutine
