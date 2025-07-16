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
subroutine pacou0(x, fvec, qt, r, c, &
                  d, fvcold, g, p, s, &
                  t, w, xold, work, check, &
                  vecr1, vecr2, typflu, vecr3, amor, &
                  masg, vecr4, vecr5, veci1, vg, &
                  indic, nbm, nmode, nt)
! aslint: disable=W1504
    implicit none
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/pacou1.h"
#include "asterfort/pacou2.h"
#include "asterfort/pacou3.h"
#include "asterfort/pacou4.h"
#include "asterfort/pacou5.h"
#include "asterfort/pacou7.h"
    integer(kind=8) :: nbm, nmode, nt
    real(kind=8) :: qt(nt, *), r(nt, *), x(*), fvec(*)
    real(kind=8) :: c(*), d(*), fvcold(*), g(*), p(*)
    real(kind=8) :: s(*), t(*), w(*), xold(*), work(*)
    real(kind=8) :: masg(*), amor(*)
    real(kind=8) :: vecr1(*), vecr2(*), vecr3(*), vecr4(*), vecr5(*)
    integer(kind=8) :: veci1(*)
    aster_logical :: restrt, sing, skip, check
    character(len=8) :: typflu
!-----------------------------------------------------------------------
    integer(kind=8) :: i, indic, its, j, k, maxits, n
    real(kind=8) :: den, eps, f, fold, stpmax, stpmx, sum
    real(kind=8) :: temp, test, tolf, tolmin, tolx, vg
!-----------------------------------------------------------------------
    parameter(eps=1.0d-8, tolx=eps, tolmin=10.d0*eps)
    parameter(tolf=1.d-04)
    parameter(maxits=200, stpmx=100.d0)
!
! FONCTION FMIN
! -------------
!
! -------------------------------------------------------------------
!
! --- TEST SEVERE POUR VOIR SI ON N'EST PAS DEJA SUR UN ZERO.
!
    check = .false.
    n = nt
!
    f = pacou2( &
        x, fvec, vecr1, vecr2, typflu, vecr3, amor, masg, vecr4, vecr5, veci1, vg, indic, nbm, &
        nmode, nt &
        )
    test = 0.0d0
    do i = 1, n
        if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
    end do
    if (test .lt. 0.01d0*tolf) goto 999
!
    sum = 0.0d0
    do i = 1, n
        sum = sum+x(i)**2
    end do
    stpmax = stpmx*max(sqrt(sum), dble(n))
    restrt = .true.
!
! --- BOUCLE PRINCIPALE.
!
    do its = 1, maxits
        if (restrt) then
!
            call pacou1(x, fvec, r, work, sqrt(eps), &
                        vecr1, vecr2, typflu, vecr3, amor, &
                        masg, vecr4, vecr5, veci1, vg, &
                        indic, nbm, nmode, nt)
            call pacou4(r, n, c, d, sing)
            if (sing) then
                check = .true.
                goto 999
            end if
            do i = 1, n
                do j = 1, n
                    qt(i, j) = 0.0d0
                end do
                qt(i, i) = 1.0d0
            end do
            do k = 1, n-1
                if (abs(c(k)) .gt. 1.0d-30) then
                    do j = 1, n
                        sum = 0.0d0
                        do i = k, n
                            sum = sum+r(i, k)*qt(i, j)
                        end do
                        sum = sum/c(k)
                        do i = k, n
                            qt(i, j) = qt(i, j)-sum*r(i, k)
                        end do
                    end do
                end if
            end do
            do i = 1, n
                r(i, i) = d(i)
                do j = 1, i-1
                    r(i, j) = 0.0d0
                end do
            end do
        else
            do i = 1, n
                s(i) = x(i)-xold(i)
            end do
            do i = 1, n
                sum = 0.0d0
                do j = 1, n
                    sum = sum+r(i, j)*s(j)
                end do
                t(i) = sum
            end do
            skip = .true.
            do i = 1, n
                sum = 0.0d0
                do j = 1, n
                    sum = sum+qt(j, i)*t(j)
                end do
                w(i) = fvec(i)-fvcold(i)-sum
                if (abs(w(i)) .ge. eps*(abs(fvec(i))+abs(fvcold(i)))) then
                    skip = .false.
!
                else
                    w(i) = 0.0d0
                end if
            end do
            if (.not. skip) then
                do i = 1, n
                    sum = 0.0d0
                    do j = 1, n
                        sum = sum+qt(i, j)*w(j)
                    end do
                    t(i) = sum
                end do
                den = 0.0d0
                do i = 1, n
                    den = den+s(i)**2
                end do
                do i = 1, n
                    s(i) = s(i)/den
                end do
!
                call pacou5(r, qt, n, t, s)
                do i = 1, n
                    if (abs(r(i, i)) .le. 1.0d-30) then
                        check = .true.
                        goto 999
                    end if
                    d(i) = r(i, i)
                end do
            end if
        end if
!
        do i = 1, n
            sum = 0.0d0
            do j = 1, n
                sum = sum+qt(i, j)*fvec(j)
            end do
            g(i) = sum
        end do
        do i = n, 1, -1
            sum = 0.0d0
            do j = 1, i
                sum = sum+r(j, i)*g(j)
            end do
            g(i) = sum
        end do
        do i = 1, n
            xold(i) = x(i)
            fvcold(i) = fvec(i)
        end do
        fold = f
        do i = 1, n
            sum = 0.0d0
            do j = 1, n
                sum = sum+qt(i, j)*fvec(j)
            end do
            p(i) = -sum
        end do
!
        call pacou7(r, n, d, p)
!
        call pacou3(xold, fold, g, p, x, &
                    f, fvec, stpmax, check, tolx, &
                    vecr1, vecr2, typflu, vecr3, amor, &
                    masg, vecr4, vecr5, veci1, vg, &
                    indic, nbm, nmode, nt)
        test = 0.0d0
        do i = 1, n
            if (abs(fvec(i)) .gt. test) test = abs(fvec(i))
        end do
        if (test .lt. tolf) then
            check = .false.
            goto 999
        end if
        if (check) then
            if (restrt) then
                goto 999
            else
                test = 0.00d0
                den = max(f, .50d0*dble(n))
                do i = 1, n
                    temp = abs(g(i))*max(abs(x(i)), 1.0d0)/den
                    if (temp .gt. test) test = temp
                end do
                if (test .lt. tolmin) then
                    goto 999
                else
                    restrt = .true.
                end if
            end if
        else
            restrt = .false.
            test = 0.0d0
            do i = 1, n
                temp = (abs(x(i)-xold(i)))/max(abs(x(i)), 1.0d0)
                if (temp .gt. test) test = temp
            end do
            if (test .lt. tolx) goto 999
        end if
    end do
    check = .true.
!
999 continue
end subroutine
