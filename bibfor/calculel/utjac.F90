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
subroutine utjac(l2d, geom, ipg, idfde, niv, &
                 ifm, nno, jacob)
! person_in_charge: olivier.boiteau at edf.fr
!-----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL LE JACOBIEN D'UN ELEMENT FINI K
!                          POUR AERER TE0003
!
! IN L2D    : FLAG INDICATEUR DU 2D
! IN GEOM   : LA GEOMETRIE
! IN IDFDE/DK/DN  : ADRESSE JEVEUX DES DERIVEES DES FONCTIONS DE FORME
! IN NIV    : NIVEAU D'IMPRESSION
! IN IFM    : UNITE LOGIQUE D'IMPRESSION
! IN NNO    : NOMBRE DE NOEUDS
! OUT JACOB : SIGNE DU JACOBIEN
!   -------------------------------------------------------------------
!     FONCTIONS INTRINSEQUES:
!       SIGN.
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       18/09/01 (OB): CREATION POUR SIMPLIFIER TE0003.F.
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
    integer(kind=8) :: ipg, idfde, niv, ifm, nno
    real(kind=8) :: jacob, geom(*)
    aster_logical :: l2d
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: i, i1, j, kp, idfdk, idfdn
    real(kind=8) :: dxde, dxdk, dyde, dydk, xp, yp, dfrde, dfrdk, dfrdn, g(3, 3)
    real(kind=8) :: j11, j21, j31
!
! INIT
    if (l2d) then
        idfdk = idfde+1
    else
        idfdn = idfde+1
        idfdk = idfdn+1
    end if
!
    if (l2d) then
! CAS 2D
        kp = 2*(ipg-1)*nno
        dxde = 0.d0
        dxdk = 0.d0
        dyde = 0.d0
        dydk = 0.d0
        do i = 1, nno
            i1 = 2*(i-1)+1
            xp = geom(i1)
            yp = geom(i1+1)
            dfrde = zr(idfde+kp+i1-1)
            dfrdk = zr(idfdk+kp+i1-1)
            dxde = dxde+xp*dfrde
            dxdk = dxdk+xp*dfrdk
            dyde = dyde+yp*dfrde
            dydk = dydk+yp*dfrdk
        end do
        jacob = dxde*dydk-dxdk*dyde
!
    else
! CAS 3D
!
        kp = 3*(ipg-1)*nno
        g(:, :) = 0.d0
        do i = 1, nno
            i1 = 3*(i-1)
            dfrde = zr(idfde+kp+i1)
            dfrdk = zr(idfdk+kp+i1)
            dfrdn = zr(idfdn+kp+i1)
            do j = 1, 3
                xp = geom(i1+j)
                g(1, j) = g(1, j)+xp*dfrde
                g(2, j) = g(2, j)+xp*dfrdn
                g(3, j) = g(3, j)+xp*dfrdk
            end do
        end do
        j11 = g(2, 2)*g(3, 3)-g(2, 3)*g(3, 2)
        j21 = g(3, 1)*g(2, 3)-g(2, 1)*g(3, 3)
        j31 = g(2, 1)*g(3, 2)-g(3, 1)*g(2, 2)
        jacob = g(1, 1)*j11+g(1, 2)*j21+g(1, 3)*j31
!
    end if
!
! CALCUL DU SIGNE DU JACOBIEN + AFFICHAGE SI NECESSAIRE
    jacob = sign(1.d0, jacob)
    if (niv .eq. 2) write (ifm, *) 'ORIENTATION MAILLE ', jacob
!
end subroutine
