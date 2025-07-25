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
subroutine proqua(quat1, quat2, quat3)
!
!
    implicit none
#include "asterfort/provec.h"
#include "blas/ddot.h"
    real(kind=8) :: quat1(4), quat2(4), quat3(4)
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (QUATERNION)
!
! CALCULE LE PRODUIT DES QUATERNIONS QUAT1 PAR QUAT2 ET MET LE
!  RESULTAT DANS QUAT3
!
! ----------------------------------------------------------------------
!
!
! IN  QUAT1  : PREMIER QUATERNION
! IN  QUAT2  : SECOND QUATERNION
! OUT QUAT3  : RESULTAT DU PRODUIT DES DEUX QUATRERNIONS
!
!
! ------------------------------------------------------------------
!
    real(kind=8) :: prosca
    integer(kind=8) :: i
    blas_int :: b_incx, b_incy, b_n
!
! ------------------------------------------------------------------
!
    call provec(quat1, quat2, quat3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    prosca = ddot(b_n, quat1, b_incx, quat2, b_incy)
    quat3(4) = quat1(4)*quat2(4)-prosca
    do i = 1, 3
        quat3(i) = quat3(i)+quat1(4)*quat2(i)+quat2(4)*quat1(i)
    end do
!
end subroutine
