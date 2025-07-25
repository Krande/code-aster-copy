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
subroutine sigeob(eps, bt, endo, ndim, lambda, &
                  mu, sigm)
!
    implicit none
#include "asterfort/diago3.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: ndim
    real(kind=8) :: eps(6), bt(6), lambda, mu
    real(kind=8) :: sigm(6), endo
!
!
!
!
!----------------------------------------------------------------------
!     CALCUL DE LA CONTRAINTE POUR LA LOI DE COMPORTEMENT
!     ENDO_ORTH_BETON
!
!     IN  EPS      : DEFORMATION
!     IN  B        : TENSEUR D ENDOMMAGEMENT DE TRACTION
!     IN  ENDO     : ENDOMMAGEMENT SCALAIRE DE COMPRESSION
!     IN  NDIM     : DIMENSION 3(3D) OU 2(2D)
!     IN  LAMBDA MU: COEFFICIENT DE LAME
!     OUT SIGM     : CONTRAINTE
!----------------------------------------------------------------------
!
    real(kind=8) :: rac2, deux, un
    real(kind=8) :: treb, treps, be(6), beeb(6), b(6)
    real(kind=8) :: to(6), tu(6), vp(3), vpe(3)
    real(kind=8) :: valbe(3), vecbe(3, 3)
    real(kind=8) :: valeps(3), veceps(3, 3), phid
    integer(kind=8) :: i, j, k, t(3, 3)
!
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
    deux = 2.d0
    rac2 = sqrt(deux)
    deux = 2.d0
    un = 1.d0
!
    phid = (un-endo)**deux
!
!      CALL DIAGO3(BT,VECB,VALB)
!      CALL R8INIR(3,0.D0,VB,1)
!      DO 32 I=1,NDIM
!          VB(I)=VALB(I)
! 32   CONTINUE
!
    call r8inir(6, 0.d0, b, 1)
!      DO 33 I=1,NDIM
!        DO 34 J=I,NDIM
!          DO 35 K=1,NDIM
!            B(T(I,J))=B(T(I,J))+VECB(I,K)*VB(K)*VECB(J,K)
! 35       CONTINUE
! 34     CONTINUE
! 33   CONTINUE
    do i = 1, 6
        b(i) = bt(i)
    end do
    call r8inir(6, 0.d0, sigm, 1)
    call r8inir(6, 0.d0, be, 1)
    do i = 1, ndim
        do j = i, ndim
            do k = 1, ndim
                be(t(i, j)) = be(t(i, j))+b(t(i, k))*eps(t(k, j))
            end do
        end do
    end do
!
    treb = 0.d0
    do i = 1, ndim
        treb = treb+be(i)
    end do
!
    treps = 0.d0
    do i = 1, ndim
        treps = treps+eps(t(i, i))
    end do
    if (treb .ge. 0.d0) then
        do i = 1, ndim
            do j = i, ndim
                sigm(t(i, j)) = sigm(t(i, j))+lambda*treb*b(t(i, j))
            end do
        end do
    end if
    if (treps .lt. 0.d0) then
        do i = 1, ndim
            sigm(t(i, i)) = sigm(t(i, i))+phid*lambda*treps
        end do
    end if
    call r8inir(6, 0.d0, beeb, 1)
    do i = 1, ndim
        do j = i, ndim
            do k = 1, ndim
                beeb(t(i, j)) = beeb(t(i, j))+b(t(i, k))*eps(t(k, j))+b(t( &
                                                                        j, k))*eps(t(k, i))
            end do
        end do
    end do
!
    call diago3(beeb, vecbe, valbe)
    call r8inir(3, 0.d0, vp, 1)
    do i = 1, ndim
        if (valbe(i) .gt. 0.d0) then
            vp(i) = valbe(i)
        else
            vp(i) = 0.d0
        end if
    end do
!
    call r8inir(6, 0.d0, to, 1)
    do i = 1, ndim
        do j = i, ndim
            do k = 1, ndim
                to(t(i, j)) = to(t(i, j))+vecbe(i, k)*vp(k)*vecbe(j, k)
            end do
        end do
    end do
!
    do i = 1, ndim
        do j = i, ndim
            do k = 1, ndim
                sigm(t(i, j)) = sigm(t(i, j))+mu/2*(to(t(i, k))*b(t(k, j))+ &
                                                    to(t(j, k))*b(t(k, i)))
            end do
        end do
    end do
    call diago3(eps, veceps, valeps)
    call r8inir(3, 0.d0, vpe, 1)
!
    do i = 1, ndim
        if (valeps(i) .lt. 0.d0) then
            vpe(i) = valeps(i)
        else
            vpe(i) = 0.d0
        end if
    end do
!
    call r8inir(6, 0.d0, tu, 1)
    do i = 1, ndim
        do j = i, ndim
            do k = 1, ndim
                tu(t(i, j)) = tu(t(i, j))+veceps(i, k)*vpe(k)*veceps(j, k)
            end do
        end do
    end do
!
    do i = 1, ndim
        do j = i, ndim
            sigm(t(i, j)) = sigm(t(i, j))+deux*mu*phid*tu(t(i, j))
        end do
    end do
!
!
    sigm(4) = rac2*sigm(4)
    sigm(5) = rac2*sigm(5)
    sigm(6) = rac2*sigm(6)
!
!
end subroutine
