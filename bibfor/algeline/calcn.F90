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
subroutine calcn(s, b, vecn)
!
    implicit none
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "blas/ddot.h"
    real(kind=8) :: b, s(6), vecn(6)
! --- BUT : CALCUL DE N ------------------------------------------------
! ======================================================================
! IN  : NDT    : NOMBRE TOTAL DE COMPOSANTES DU TENSEUR ----------------
! --- : NDI    : NOMBRE DE COMPOSANTES DIAGONALES DU TENSEUR -----------
! --- : S      : DEVIATEUR DES CONTRAINTES -----------------------------
! --- : B      : PARAMETRE DU CALCUL DE LA NORMALE ---------------------
! OUT : VECN   : N = (B*S/SII+I)/SQRT(B**2+3) --------------------------
! ======================================================================
    integer(kind=8) :: ii, ndt, ndi
    real(kind=8) :: sii, racine, un, trois
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
! --- INITIALISATION DE PARAMETRE --------------------------------------
! ======================================================================
    parameter(un=1.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
    call jemarq()
! ======================================================================
! --- INITIALISATION ---------------------------------------------------
! ======================================================================
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    sii = ddot(b_n, s, b_incx, s, b_incy)
    sii = sqrt(sii)
! ======================================================================
! --- CALCUL DE N ------------------------------------------------------
! ======================================================================
    racine = sqrt(b*b+trois)
    do ii = 1, ndt
        vecn(ii) = b*s(ii)/(sii*racine)
    end do
    do ii = 1, ndi
        vecn(ii) = vecn(ii)+un/racine
    end do
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
