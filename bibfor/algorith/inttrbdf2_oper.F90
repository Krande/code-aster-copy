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

subroutine inttrbdf2_oper(nbequ, par, mgen, kgen, cgen, &
                          ktilda, ftild1, ftild2)
    implicit none
!
! person_in_charge: nicolas.tardieu at edf.fr
!
! inttrbdf2_oper : Calculate (or update) the operators for TR-BDF2 integration
!
!             --- ktilda, ftild1, ftild2 ---
!           kt (i,j) = a1*m(i,j) + k(i,j) + a2*c(i,j)
!           ft1(i,j) = a2*m(i,j)
!           ft2(i,j) = a1*m(i,j) + a2*c(i,j)
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/trlds.h"
!
!   -0.1- Input/output arguments
    integer(kind=8), intent(in)           :: nbequ
    real(kind=8)                       :: par(:)
    real(kind=8), pointer  :: mgen(:), kgen(:), cgen(:)
    real(kind=8), pointer :: ktilda(:), ftild1(:), ftild2(:)
!
!   -0.2- Local variables
    integer(kind=8)           :: i, j, iret

!   --------------------------------------------------------------------------
#define a(k) par(k)
#define norm_coef par(9)
#define mdiag_r par(12)
#define kdiag_r par(13)
#define cdiag_r par(14)

#define mdiag (nint(mdiag_r).eq.1)
#define kdiag (nint(kdiag_r).eq.1)
#define cdiag (nint(cdiag_r).eq.1)

#define m(row,col) mgen((col-1)*nbequ+row)
#define k(row,col) kgen((col-1)*nbequ+row)
#define c(row,col) cgen((col-1)*nbequ+row)

#define kt(row,col) ktilda((col-1)*nbequ+row)
#define ft1(row,col) ftild1((col-1)*nbequ+row)
#define ft2(row,col) ftild2((col-1)*nbequ+row)

!       --- M is diagonal
    if (mdiag) then
        if (kdiag) then
            if (cdiag) then
                do i = 1, nbequ
                    ftild2(i) = a(1)*mgen(i)+a(2)*cgen(i)
                    ktilda(i) = ftild2(i)+kgen(i)
                    ftild1(i) = a(2)*mgen(i)
                end do
            else
                do i = 1, nbequ
                    ft2(i, i) = a(1)*mgen(i)+a(2)*c(i, i)
                    kt(i, i) = ft2(i, i)+kgen(i)
                    ft1(i, i) = a(2)*mgen(i)
                    do j = i+1, nbequ
                        kt(i, j) = a(2)*c(i, j)
                        kt(j, i) = a(2)*c(j, i)
                        ft2(i, j) = kt(i, j)
                        ft2(j, i) = kt(j, i)
                    end do
                end do
            end if
        else
            if (cdiag) then
                do i = 1, nbequ
                    ftild2(i) = a(1)*mgen(i)+a(2)*cgen(i)
                    kt(i, i) = ftild2(i)+k(i, i)
                    ftild1(i) = a(2)*mgen(i)
                    do j = i+1, nbequ
                        kt(i, j) = k(i, j)
                        kt(j, i) = k(j, i)
                    end do
                end do
            else
                do i = 1, nbequ
                    ft2(i, i) = a(1)*mgen(i)+a(2)*c(i, i)
                    kt(i, i) = ft2(i, i)+k(i, i)
                    ft1(i, i) = a(2)*mgen(i)
                    do j = i+1, nbequ
                        ft2(i, j) = a(2)*c(i, j)
                        ft2(j, i) = a(2)*c(j, i)
                        kt(i, j) = ft2(i, j)+k(i, j)
                        kt(j, i) = ft2(j, i)+k(j, i)
                    end do
                end do
            end if
        end if

!   --- M is not diagonal, K is supposed to be full as well
    else
        if (cdiag) then
            do i = 1, nbequ
                ft2(i, i) = a(1)*m(i, i)+a(2)*cgen(i)
                kt(i, i) = ft2(i, i)+k(i, i)
                ft1(i, i) = a(2)*m(i, i)
                do j = i+1, nbequ
                    ft2(i, j) = a(1)*m(i, j)
                    ft2(j, i) = a(1)*m(j, i)
                    kt(i, j) = ft2(i, j)+k(i, j)
                    kt(j, i) = ft2(j, i)+k(j, i)
                    ft1(i, j) = a(2)*m(i, j)
                    ft1(j, i) = a(2)*m(j, i)
                end do
            end do
        else
            do i = 1, nbequ
                ft2(i, i) = a(1)*m(i, i)+a(2)*c(i, i)
                kt(i, i) = ft2(i, i)+k(i, i)
                ft1(i, i) = a(2)*m(i, i)
                do j = i+1, nbequ
                    ft2(i, j) = a(1)*m(i, j)+a(2)*c(i, j)
                    ft2(j, i) = a(1)*m(j, i)+a(2)*c(j, i)
                    kt(i, j) = ft2(i, j)+k(i, j)
                    kt(j, i) = ft2(j, i)+k(j, i)
                    ft1(i, j) = a(2)*m(i, j)
                    ft1(j, i) = a(2)*m(j, i)
                end do
            end do
        end if
    end if

    norm_coef = -1.d25
    do i = 1, size(ktilda)
        if (abs(ktilda(i)) .gt. norm_coef) norm_coef = abs(ktilda(i))
    end do

    ASSERT(norm_coef .gt. 1.d-25)
    do i = 1, size(ftild1)
        ftild1(i) = ftild1(i)/norm_coef
        ftild2(i) = ftild2(i)/norm_coef
    end do

    do i = 1, size(ktilda)
        ktilda(i) = ktilda(i)/norm_coef
    end do

!   --- Factorize ktilda if needed for later resolution for displacement
    if (size(ktilda) .gt. nbequ) then
        call trlds(ktilda, nbequ, nbequ, iret)
    end if

end subroutine
