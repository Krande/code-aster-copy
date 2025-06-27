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
subroutine postsm(option, fm, df, sigm, sigp, &
                  dsidep)
!
!
    implicit none
#include "asterf_types.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    character(len=16) :: option
    real(kind=8) :: fm(3, 3), df(3, 3), sigm(6), sigp(6), dsidep(6, 3, 3)
!
!------------------------------------------------------------
!   IN  OPTION : OPTION DEMANDEE : RIGI_MECA_TANG, FULL_MECA, RAPH_MECA
!   IN  FM : GRADIENT DE LA TRANSFORMATION EN T-
!   IN  DF : GRADIENT DE LA TRANSFORMATION DE T- A T+
!   IN  SIGM : CONTRAINTE EN T-
!   IN/OUT  SIGP : CONTRAINTE CAUCHY EN T+ -> CONTRAINTE KIRCHHOF EN T+
!   IN/OUT  DSIDEP : MATRICE TANGENTE D(SIG)/DF  ->
!                    D(TAU)/D(FD) * (FD)t
!-----------------------------------------------------------------------
    aster_logical :: resi, rigi
    integer(kind=8) :: kl, p, q, i
    real(kind=8) :: jm, dj, jp, tau(6), j, mat(6, 3, 3), id(3, 3), rc(6)
!
    data id/1.d0, 0.d0, 0.d0,&
     &              0.d0, 1.d0, 0.d0,&
     &              0.d0, 0.d0, 1.d0/
!
    real(kind=8) :: rac2
    blas_int :: b_incx, b_incy, b_n
    parameter(rac2=sqrt(2.d0))
    data rc/1.d0, 1.d0, 1.d0, rac2, rac2, rac2/
!
!
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
!
    jm = fm(1, 1)*(fm(2, 2)*fm(3, 3)-fm(2, 3)*fm(3, 2))-fm(2, 1)*(fm(1, 2)*fm(3, 3)-fm(1, 3)*fm(3&
         &, 2))+fm(3, 1)*(fm(1, 2)*fm(2, 3)-fm(1, 3)*fm(2, 2))
!
    dj = df(1, 1)*(df(2, 2)*df(3, 3)-df(2, 3)*df(3, 2))-df(2, 1)*(df(1, 2)*df(3, 3)-df(1, 3)*df(3&
         &, 2))+df(3, 1)*(df(1, 2)*df(2, 3)-df(1, 3)*df(2, 2))
!
    jp = jm*dj
!
    if (resi) then
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        call dscal(b_n, jp, sigp, b_incx)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigp, b_incx, tau, b_incy)
        j = jp
    else
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigm, b_incx, tau, b_incy)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        call dscal(b_n, jm, tau, b_incx)
        j = jm
    end if
!
!
    if (rigi) then
        b_n = to_blas_int(54)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dsidep, b_incx, mat, b_incy)
        b_n = to_blas_int(54)
        b_incx = to_blas_int(1)
        call dscal(b_n, j, mat, b_incx)
        do kl = 1, 6
            do p = 1, 3
                do q = 1, 3
                    dsidep(kl, p, q) = tau(kl)*id(p, q)
                    do i = 1, 3
                        dsidep(kl, p, q) = dsidep(kl, p, q)+mat(kl, p, i)*df(q, i)
                    end do
!
                    dsidep(kl, p, q) = dsidep(kl, p, q)*rc(kl)
!
                end do
            end do
        end do
    end if
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, rac2, sigp(4), b_incx)
!
end subroutine
