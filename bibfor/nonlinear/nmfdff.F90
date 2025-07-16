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
subroutine nmfdff(ndim, nno, axi, g, r, &
                  rigi, matsym, fr, vff, dff, &
                  def, pff)
!     BUT:  CALCUL DES PRODUITS DU GRADIENT DE TRANSFORMATION FR
!     PAR LES DERIVEES DE FONCTIONS DE FORME DFF POUR NMDLOG
!     CONFIGURATION LAGRANGIENNE.
! ----------------------------------------------------------------------
! IN  NDIM    : DIMENSION DU PROBLEME : 2 OU 3
! IN  DFF     : DERIVEE DES FONCTIONS DE FORME
! IN  FR      : GRADIENT TRANSFORMATION
! OUT DEF     : PRODUITS F*DFF
! OUT PFF     : PRODUITS DFF*DFF
!
    implicit none
#include "asterf_types.h"
    integer(kind=8) :: ndim, nno, n, nmax, i, m, g
    real(kind=8) :: dff(nno, *), vff(nno, *), fr(3, 3), rac2, r
    real(kind=8) :: def(2*ndim, nno, ndim), pff(2*ndim, nno, nno)
    aster_logical :: axi, rigi, matsym
!
    rac2 = sqrt(2.d0)
!
    if (ndim .eq. 3) then
!
        do n = 1, nno
            do i = 1, 3
                def(1, n, i) = fr(i, 1)*dff(n, 1)
                def(2, n, i) = fr(i, 2)*dff(n, 2)
                def(3, n, i) = fr(i, 3)*dff(n, 3)
                def(4, n, i) = (fr(i, 1)*dff(n, 2)+fr(i, 2)*dff(n, 1))/rac2
                def(5, n, i) = (fr(i, 1)*dff(n, 3)+fr(i, 3)*dff(n, 1))/rac2
                def(6, n, i) = (fr(i, 2)*dff(n, 3)+fr(i, 3)*dff(n, 2))/rac2
            end do
        end do
!
        if (rigi) then
            do n = 1, nno
                if (matsym) then
                    nmax = n
                else
                    nmax = nno
                end if
                do m = 1, nmax
                    pff(1, n, m) = dff(n, 1)*dff(m, 1)
                    pff(2, n, m) = dff(n, 2)*dff(m, 2)
                    pff(3, n, m) = dff(n, 3)*dff(m, 3)
                    pff(4, n, m) = (dff(n, 1)*dff(m, 2)+dff(n, 2)*dff(m, 1))/ &
                                   rac2
                    pff(5, n, m) = (dff(n, 1)*dff(m, 3)+dff(n, 3)*dff(m, 1))/ &
                                   rac2
                    pff(6, n, m) = (dff(n, 2)*dff(m, 3)+dff(n, 3)*dff(m, 2))/ &
                                   rac2
                end do
            end do
        end if
!
    else if (ndim .eq. 2) then
!
        do n = 1, nno
            do i = 1, 2
                def(1, n, i) = fr(i, 1)*dff(n, 1)
                def(2, n, i) = fr(i, 2)*dff(n, 2)
                def(3, n, i) = 0.d0
                def(4, n, i) = (fr(i, 1)*dff(n, 2)+fr(i, 2)*dff(n, 1))/rac2
            end do
        end do
! 5.2.5 - TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
        if (axi) then
            do n = 1, nno
                def(3, n, 1) = fr(3, 3)*vff(n, g)/r
            end do
        end if
        if (rigi) then
            do n = 1, nno
                if (matsym) then
                    nmax = n
                else
                    nmax = nno
                end if
                do m = 1, nmax
                    pff(1, n, m) = dff(n, 1)*dff(m, 1)
                    pff(2, n, m) = dff(n, 2)*dff(m, 2)
                    pff(3, n, m) = 0.d0
                    pff(4, n, m) = (dff(n, 1)*dff(m, 2)+dff(n, 2)*dff(m, 1))/ &
                                   rac2
                end do
            end do
        end if
!
    end if
end subroutine
