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
subroutine ordon1(vale, nb)
! aslint: disable=W1306
    implicit none
#include "asterfort/ordr8.h"
#include "blas/dcopy.h"
    integer(kind=8) :: nb
    real(kind=8) :: vale(*)
! person_in_charge: mathieu.courtois at edf.fr
! ----------------------------------------------------------------------
!     APPELEE PAR ORDONN : ON SAIT DEJA QU'IL FAUT INVERSER L'ORDRE
!     ORDON1 POUR LES FONCTIONS A VALEURS REELLES
! IN/OUT : VALE : ABSCISSES SUIVIES DES ORDONNEES
! IN     : NB   : NBRE DE POINTS
! ----------------------------------------------------------------------
    integer(kind=8) :: i, iord(nb)
    real(kind=8) :: xbid(nb), yrbid(nb)
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
!
!
    b_n = to_blas_int(nb)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vale, b_incx, xbid, b_incy)
    b_n = to_blas_int(nb)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vale(nb+1), b_incx, yrbid, b_incy)
    call ordr8(xbid, nb, iord)
    do i = 1, nb
        vale(i) = xbid(iord(i))
        vale(nb+i) = yrbid(iord(i))
    end do
!
end subroutine
