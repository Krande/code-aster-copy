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
! aslint: disable=W1306
!
subroutine reerel(elrefp, nnop, ndim, tabar, xe, &
                  xg, dxg)
!
    implicit none
!
#include "asterfort/elrfvf.h"
#include "asterfort/elrfdf.h"
!
    integer(kind=8), intent(in) :: ndim, nnop
    real(kind=8), intent(in) :: xe(ndim), tabar(*)
    character(len=8), intent(in) :: elrefp
    real(kind=8), intent(out) :: xg(ndim)
    real(kind=8), intent(out), optional :: dxg(ndim)
!
!
!         TROUVER LES COORDONNEES REELLES D'UN POINT
!         A PARTIR DE SES COORDONNEES DE REFERENCE
!
!     ENTREE
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       IGEOM   : COORDONNEES DES NOEUDS DE L'ELEMENT
!       ELP     : TYPE DE L'ELEMENT
!       XE      : COORDONNES DE REFERENCE DU POINT
!
!     SORTIE
!       XG       : COORDONNES REELLES DU POINT
!       DXG      : DERIVEES DES COORDONNES REELLES DU POINT
!......................................................................
!
    real(kind=8) :: ff(nnop), dff(nnop)
    integer(kind=8) :: i, j
!
!......................................................................
!
    xg = 0.d0
    if (present(dxg)) then
        dxg = 0.d0
    end if
!
! --- VALEURS DES FONCTIONS DE FORME EN XE: FF
!
    if (elrefp(1:2) .eq. 'SE') then
        call elrfvf(elrefp, xe(1), ff)
    else
        call elrfvf(elrefp, xe, ff)
    end if
    if (present(dxg)) then
        if (elrefp(1:2) .eq. 'SE') then
            call elrfdf(elrefp, xe(1), dff)
        else
            call elrfdf(elrefp, xe, dff)
        end if
    end if
!
! --- COORDONNES DU POINT DANS L'ELEMENT REEL
!
    do j = 1, ndim
        do i = 1, nnop
            xg(j) = xg(j)+tabar(ndim*(i-1)+j)*ff(i)
        end do
    end do
!
    if (present(dxg)) then
        do j = 1, ndim
            do i = 1, nnop
                dxg(j) = dxg(j)+tabar(ndim*(i-1)+j)*dff(i)
            end do
        end do
    end if
!
end subroutine
