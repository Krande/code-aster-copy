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
subroutine xcedge(ndime, pinref, pi1, pi2, pmiref, &
                  m12, crit)
!
    implicit none
!
#include "asterfort/xnormv.h"
#include "blas/ddot.h"
!
    real(kind=8) :: pinref(*), pmiref(*), crit
    integer(kind=8) :: pi1, pi2, m12, ndime
!
!
! FONCTION REALISEE:  CALCUL D UN CRITERE DE COURBURE SUR UNE ARETE DANS
!                          L ELEMENT DE REFERENCE
!   CRIT=0 SI ARETE DROITE
!       >0 SI ARETE COURBE
!
    real(kind=8) :: pipk(ndime), pimik(ndime), rbid, cosi
    integer(kind=8) :: i
    blas_int :: b_incx, b_incy, b_n
!
    crit = 0.d0
    do i = 1, ndime
        pipk(i) = pinref(ndime*(pi2-1)+i)-pinref(ndime*(pi1-1)+i)
        pimik(i) = pmiref(ndime*(m12-1)+i)-pinref(ndime*(pi1-1)+i)
    end do
    call xnormv(ndime, pipk, rbid)
    call xnormv(ndime, pimik, rbid)
    b_n = to_blas_int(ndime)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    cosi = ddot(b_n, pipk, b_incx, pimik, b_incy)
!    write(6,*)'xcedge: cosi=', cosi
    if (cosi .lt. 1.d0) crit = dacos(cosi)
!
end subroutine
