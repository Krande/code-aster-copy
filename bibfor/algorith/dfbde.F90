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
subroutine dfbde(dim, b, e, deuxmu, lambda, &
                 dsidep)
!
!
    implicit none
#include "asterfort/dfpdf.h"
#include "asterfort/diago3.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: dim
    real(kind=8) :: b(6), e(6), deuxmu, lambda
    real(kind=8) :: dsidep(6, 6)
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT ENDO_ORTH_BETON
!     CALCUL DE LA DERIVEE DE LA FORCE THERMODYNAMIQUE (ENDO TRACTION)
!     PAR RAPPORT A LA DEFORMATION:DFB/DEPS
!
!     FB=-LAMBDA.TR(EB).H(TR(EB))-MU/2*((BE+EB)_+*E + E*(BE+EB)_+)
!        +ECROB*(I-B)
!     IN  DIM      : DIMENSION 3(3D) OU 2(2D)
!     IN  E        : DEFORMATION
!     IN  B        : TENSEUR D ENDOMMAGEMENT DE TRACTION
!     IN  LAMBDA   : /
!     IN  DEUXMU   : / COEFFICIENTS DE LAME
!     OUT DSIDEP      : DFB/DEPS
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, l, t(3, 3)
    integer(kind=8) :: p, q, ik, pq, rs
    real(kind=8) :: rtemp2
    real(kind=8) :: rtemp3, rtemp4, rac2, kron(3, 3), dbede(6, 6), mtemp(6, 6)
    real(kind=8) :: beeb(6), beebp(6), vecbeb(3, 3), valbeb(3)
    real(kind=8) :: f1b(6, 6), f2b(6, 6)
    real(kind=8) :: a(6, 6), treb, c(6, 6), trec
!
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
!
    kron(1, 1) = 1.d0
    kron(1, 2) = 0.d0
    kron(1, 3) = 0.d0
    kron(2, 1) = 0.d0
    kron(2, 2) = 1.d0
    kron(2, 3) = 0.d0
    kron(3, 1) = 0.d0
    kron(3, 2) = 0.d0
    kron(3, 3) = 1.d0
!
    rac2 = sqrt(2.d0)
!
    call r8inir(36, 0.d0, dbede, 1)
    call r8inir(36, 0.d0, c, 1)
    call r8inir(36, 0.d0, dsidep, 1)
    call r8inir(36, 0.d0, a, 1)
    call r8inir(6, 0.d0, beeb, 1)
    call r8inir(36, 0.d0, f1b, 1)
    call r8inir(36, 0.d0, f2b, 1)
!
!-----CALCUL DE D(BE+EB)/DE------------------------------------
    rtemp3 = 0.d0
    rtemp4 = 0.d0
    do i = 1, dim
        do j = i, dim
            do k = 1, dim
                do l = 1, dim
                    if (i .eq. j) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = rac2
                    end if
                    if (k .eq. l) then
                        rtemp4 = 1.d0
                    else
                        rtemp4 = 1.d0/rac2
                    end if
                    dbede(t(i, j), t(k, l)) = dbede(t(i, j), t(k, l))+(kron(k, &
                                                  i)*b(t(l, j))+kron(j, l)*b(t(i, k)))*rtemp3*rtemp4
                end do
                beeb(t(i, j)) = beeb(t(i, j))+b(t(i, k))*e(t(k, j))+e(t(i, &
                                                                        k))*b(t(k, j))
            end do
        end do
    end do
!
!
    call diago3(beeb, vecbeb, valbeb)
    do i = 1, 3
        if (valbeb(i) .lt. 0.d0) then
            valbeb(i) = 0.d0
        end if
    end do
!
    call r8inir(6, 0.d0, beebp, 1)
!
    do i = 1, 3
        do j = i, 3
            do k = 1, 3
                beebp(t(i, j)) = beebp(t(i, j))+vecbeb(i, k)*valbeb(k)* &
                                 vecbeb(j, k)
            end do
        end do
    end do
!
!
!
!
!----------------------------------------------------------------------
!
!-----CALCUL DE F2B----------------------------------------------------
!
!
    call dfpdf(6, beeb, mtemp)
!
    do ik = 1, 6
        do pq = 1, 6
            do rs = 1, 6
                a(ik, pq) = a(ik, pq)+mtemp(ik, rs)*dbede(rs, pq)
            end do
        end do
    end do
!
    rtemp2 = 0.d0
    rtemp3 = 0.d0
    rtemp4 = 0.d0
    do i = 1, 3
        do j = i, 3
            do k = 1, 3
                do l = 1, 6
                    if (i .eq. j) then
                        rtemp2 = 1.d0
                    else
                        rtemp2 = rac2
                    end if
                    if (k .eq. j) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = 1.d0/rac2
                    end if
                    if (k .eq. i) then
                        rtemp4 = 1.d0
                    else
                        rtemp4 = 1.d0/rac2
                    end if
                    f2b(t(i, j), l) = f2b(t(i, j), l)+(a(t(i, k), l)*e(t(k, j)) &
                                                     *rtemp4+e(t(i, k))*a(t(k, j), l)*rtemp3)*rtemp2
                end do
            end do
        end do
    end do
!
!
    rtemp2 = 0.d0
    rtemp3 = 0.d0
    do i = 1, 3
        do j = i, 3
            do p = 1, 3
                do q = 1, 3
                    if (i .eq. j) then
                        rtemp2 = 1.d0
                    else
                        rtemp2 = rac2
                    end if
                    if (p .eq. q) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = 1.d0/rac2
                    end if
                    f2b(t(i, j), t(p, q)) = f2b(t(i, j), t(p, q))+(kron(j, q)* &
                                                 beebp(t(i, p))+kron(i, p)*beebp(t(j, q)))*rtemp2* &
                                            rtemp3
                end do
            end do
        end do
    end do
!
!-----------------------------------------------------------------------
!
!-----CALCUL DE F1B-----------------------------------------------------
!
    treb = (beeb(1)+beeb(2)+beeb(3))/2.d0
!
    if (treb .lt. 0) then
        trec = 0.d0
    else
        trec = 1.d0
    end if
!
!
    rtemp2 = 0.d0
    rtemp3 = 0.d0
    do i = 1, 3
        do j = i, 3
            do p = 1, 3
                do q = 1, 3
                    if (i .eq. j) then
                        rtemp2 = 1.d0
                    else
                        rtemp2 = rac2
                    end if
                    if (p .eq. q) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = 1.d0/rac2
                    end if
                    f1b(t(i, j), t(p, q)) = f1b(t(i, j), t(p, q))+trec*b(t(p, &
                                                                        q))*e(t(i, j))*rtemp2*rtemp3
                end do
            end do
        end do
    end do
!
!
    rtemp2 = 0.d0
    rtemp3 = 0.d0
    do i = 1, 3
        do j = i, 3
            do p = 1, 3
                do q = 1, 3
                    if (i .eq. j) then
                        rtemp2 = 1.d0
                    else
                        rtemp2 = rac2
                    end if
                    if (p .eq. q) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = 1.d0/rac2
                    end if
                    f1b(t(i, j), t(p, q)) = f1b(t(i, j), t(p, q))+trec*treb &
                                            *kron(i, p)*kron(j, q)*rtemp2*rtemp3
                end do
            end do
        end do
    end do
!
!
!----CALCUL DE FB-------------------------------------------------------
!
    do i = 1, 6
        do j = 1, 6
            dsidep(i, j) = -lambda*f1b(i, j)-deuxmu/4.d0*f2b(i, j)
        end do
    end do
!
end subroutine
