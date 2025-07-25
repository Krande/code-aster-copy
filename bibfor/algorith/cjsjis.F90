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

subroutine cjsjis(mod, mater, deps, yd, yf, &
                  r, drdy)
    implicit none
!     INTEGRATION PLASTIQUE (MECANISME ISOTROPE SEUL) DE LA LOI CJS
!
!     RESOLUTION PAR METHODE DE NEWTON       DRDY(DYI) DDYI = - R(DYI)
!
!     CALCUL DU SECOND MEMBRE : - R(DYI)
!     CALCUL DU JACOBIEN      : DRDY(DYI)
!                               Y =  ( SIG    ,  VIN    ,  LAMBI   )
!                               R = -( LE     ,  LQ     ,  FI      )
!                           DRDY  =  ( DLEDS  ,  DLEDQ  ,  DLEDL   )
!                                    ( DLQDS  ,  DLQDQ  ,  DLQDL   )
!                                    ( DFIDS  ,  DFIDQ  ,  DFIDL   )
!     ------------------------------------------------------------------
!     IN   MOD      :  MODELISATION
!          MATER    :  COEFFICIENTS MATERIAU A T+DT
!          DEPS     :  INCREMENT DE DEFORMATION
!          YD       :  VARIABLES A T = (SIGD, VIND, LAMBID)
!          YF       :  VARIABLES A T+DT = (SIGF, VINF, LAMBIF)
!     OUT  R        :  SECOND MEMBRE
!          DRDY     :  JACOBIEN
!     ------------------------------------------------------------------
!
#include "asterfort/lcicma.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nmod
    parameter(nmod=8)
!
    real(kind=8) :: deps(6)
    real(kind=8) :: yd(nmod), yf(nmod), r(nmod), drdy(nmod, nmod)
    real(kind=8) :: mater(14, 2), n, kop, pa
    real(kind=8) :: hooknl(6, 6), hook(6, 6)
    real(kind=8) :: e, nu, al, la, mu, i1f
    real(kind=8) :: dlambi, coef1, coef2
    real(kind=8) :: le(6), lq, fi
    real(kind=8) :: dleds(6, 6), dledq(6), dledl(6)
    real(kind=8) :: dlqds(6), dlqdq, dlqdl
    real(kind=8) :: dfids(6), dfidq, dfidl
    real(kind=8) :: dsignl(6), dsigl(6), depse(6), qinit
    integer(kind=8) :: i, j
!
    real(kind=8) :: zero, un, d12, deux, trois, kron(6), iden6(6, 6)
!
    parameter(d12=.5d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
!
!
    character(len=8) :: mod
!
    common/tdim/ndt, ndi
!
    data iden6/un, zero, zero, zero, zero, zero,&
         &                   zero, un, zero, zero, zero, zero,&
         &                   zero, zero, un, zero, zero, zero,&
         &                   zero, zero, zero, un, zero, zero,&
         &                   zero, zero, zero, zero, un, zero,&
         &                   zero, zero, zero, zero, zero, un/
!
    data kron/un, un, un, zero, zero, zero/
!
!
!-----------------------------------------------------------------------
!->     PROPRIETES CJS MATERIAU
!------------------------------
!
    n = mater(3, 2)
    kop = mater(4, 2)
    pa = mater(12, 2)
    qinit = mater(13, 2)
!
!-----------------------------------------------------------------------
!->     OPERATEURS DE RIGIDITE
!-----------------------------------------------------------------------
!
!
!- OPERATEUR LINEAIRE
!++++++++++++++++++++
!
    hook(:, :) = zero
!
    e = mater(1, 1)
    nu = mater(2, 1)
    al = e*(un-nu)/(un+nu)/(un-deux*nu)
    la = nu*e/(un+nu)/(un-deux*nu)
    mu = e*d12/(un+nu)
!
!
! - 3D/DP/AX
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        do i = 1, ndi
            do j = 1, ndi
                if (i .eq. j) hook(i, j) = al
                if (i .ne. j) hook(i, j) = la
            end do
        end do
        do i = ndi+1, ndt
            do j = ndi+1, ndt
                if (i .eq. j) hook(i, j) = deux*mu
            end do
        end do
        !
        ! - CP/1D
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'ALGORITH2_15')
    end if
!
!
!- OPERATEUR NON LINEAIRE
!++++++++++++++++++++++++
!
    i1f = zero
    do i = 1, ndi
        i1f = i1f+yf(i)
    end do
    if ((i1f+qinit) .eq. 0.d0) then
        i1f = -qinit+1.d-12*pa
    end if
!
    coef1 = ((i1f+qinit)/trois/pa)**n
!
    do i = 1, ndt
        do j = 1, ndt
            hooknl(i, j) = coef1*hook(i, j)
        end do
    end do
!
!
!
!
!
!
!--------------------------------------------------------
!->     LOI D ETAT : LE
!->  ET DERIVEE DE LA LOI D ETAT : DLEDS, DLEDQ, DLEDL
!--------------------------------------------------------
!
!- LOI D ETAT
!++++++++++++
!
    dlambi = yf(ndt+2)-yd(ndt+2)
!
    do i = 1, ndt
        depse(i) = deps(i)+dlambi*kron(i)/trois
    end do
!
    dsignl(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), depse(1:ndt))
!
    do i = 1, ndt
        le(i) = yf(i)-yd(i)-dsignl(i)
    end do
!
!
!- DERIVEE DE LA LOI D ETAT
!++++++++++++++++++++++++++
!
    coef2 = n/trois/pa*(trois*pa/(i1f+qinit))**(un-n)
    dsigl(1:ndt) = matmul(hook(1:ndt, 1:ndt), depse(1:ndt))
    dleds(:, :) = zero
!
    do i = 1, ndt
        !
        do j = 1, ndt
            dleds(i, j) = iden6(i, j)-coef2*dsigl(i)*kron(j)
        end do
        !
        dledq(i) = zero
        !
        dledl(i) = zero
        do j = 1, ndt
            dledl(i) = dledl(i)-hooknl(i, j)*kron(j)/trois
        end do
        !
    end do
!
!
!
!
!
!-----------------------------------------------------------------------
!->     LOI D ECROUISSAGE DE QISO : LQ
!->  ET DERIVEE DE LA LOI D ECROUISSAGE DE QISO : DLQDS, DLQDQ, DLQDL
!-----------------------------------------------------------------------
!
!
!- LOI D ECROUISSAGE
!+++++++++++++++++++
    lq = yf(ndt+1)-yd(ndt+1)+dlambi*kop*(yf(ndt+1)/pa)**n
!
!
!- DERIVEE DE LA LOI D ECROUISSAGE
!+++++++++++++++++++++++++++++++++
!
    dlqds(:) = zero
    dlqdq = un+dlambi*kop*n/pa*(pa/yf(ndt+1))**(un-n)
    dlqdl = kop*(yf(ndt+1)/pa)**n
!------------------------------------------------------------------
!->     SEUIL ISOTROPE : FI
!->  ET DERIVEE DE LA FONCTION SEUIL ISOTROPE : DFIDS, DFIDQ, DFIDL
!------------------------------------------------------------------
!
!
!- SEUIL ISOTROPE
!++++++++++++++++
!
    fi = -(i1f+qinit)/trois+yf(ndt+1)
!
!
!- DERIVEE DU SEUIL ISOTROPE
!+++++++++++++++++++++++++++
!
    do i = 1, ndt
        dfids(i) = -kron(i)/trois
    end do
!
    dfidq = un
    dfidl = zero
!
!
!
!
!-------------------------------------------
!->     ASSEMBLAGE DE R = - ( LE, LQ, FI )
!->  ET ASSEMBLAGE DE DRDY
!
!       DRDY  =   DLEDS   DLEDQ   DLEDL
!                 DLQDS   DLQDQ   DLQDL
!                 DFIDS   DFIDQ   DFIDL
!
!-------------------------------------------
!
!
!- ASSEMBLAGE DE R
!+++++++++++++++++
!
    do i = 1, ndt
        r(i) = -le(i)
    end do
!
    r(ndt+1) = -lq
    r(ndt+2) = -fi
!
!
!- ASSEMBLAGE DE DRDY
!++++++++++++++++++++
!
    call lcicma(dleds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                1, 1)
    call lcicma(dledq, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+1)
    call lcicma(dledl, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+2)
!
    call lcicma(dlqds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+1, 1)
    drdy(ndt+1, ndt+1) = dlqdq
    drdy(ndt+1, ndt+2) = dlqdl
!
    call lcicma(dfids, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+2, 1)
    drdy(ndt+2, ndt+1) = dfidq
    drdy(ndt+2, ndt+2) = dfidl
!
!
end subroutine
