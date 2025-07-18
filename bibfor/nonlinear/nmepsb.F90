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
subroutine nmepsb(ndim, nno, axi, vff, dfdi, &
                  deplg, epsb, geps)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    aster_logical :: axi
    integer(kind=8) :: ndim, nno
    real(kind=8) :: vff(nno), dfdi(nno, ndim), deplg(*)
    real(kind=8) :: epsb(6), geps(6, 3)
!
! ----------------------------------------------------------------------
!       CALCUL DES DEFORMATIONS REGULARISEES ET LEURS GRADIENTS
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO     : NOMBRE DE NOEUDS (FAMILLE E-BARRE)
! IN  AXI     : .TRUE. SI AXISYMETRIQUE
! IN  VFF     : VALEURS DES FONCTIONS DE FORME (POINT DE GAUSS COURANT)
! IN  DFDI    : DERIVEE DES FONCTIONS DE FORME (POINT DE GAUSS COURANT)
! IN  DEPLG   : DEPLACEMENT GENERALISE (U ET E-BARRE)
! OUT EPSB    : DEFORMATIONS REGULARISEES EPSB(6)
! OUT GEPS    : GRADIENT DES DEFORMATIONS REGULARISEES GEPS(6,3)
! ----------------------------------------------------------------------
!
    integer(kind=8) :: kl, i, ndimsi, ndl
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
!
!
!
    ndimsi = 2*ndim
    ndl = ndim+ndimsi
!
    call r8inir(6, 0.d0, epsb, 1)
    call r8inir(18, 0.d0, geps, 1)
    do kl = 1, ndimsi
        b_n = to_blas_int(nno)
        b_incx = to_blas_int(ndl)
        b_incy = to_blas_int(1)
        epsb(kl) = ddot(b_n, deplg(kl+ndim), b_incx, vff, b_incy)
        do i = 1, ndim
            b_n = to_blas_int(nno)
            b_incx = to_blas_int(ndl)
            b_incy = to_blas_int(1)
            geps(kl, i) = ddot(b_n, deplg(kl+ndim), b_incx, dfdi(1, i), b_incy)
        end do
    end do
!
    if (axi) then
        call utmess('F', 'ALGORITH7_76')
    end if
end subroutine
