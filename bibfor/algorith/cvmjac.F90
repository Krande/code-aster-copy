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

subroutine cvmjac(mod, nmat, materf, timed, timef, &
                  yf, dy, nmod, epsd, deps, &
                  drdy)
! aslint: disable=W1501
    implicit none
!       VISCOCHABOCHE   :
!           CALCUL DU JACOBIEN DU SYSTEME NL A RESOUDRE = DRDY
!                   DY  = ( DSIG DX1 DX2 DP DR 0 (DEPS3))
!                   Y   = ( SIG  X1  X2  P  R  Q (EPS3) )
!
!       DRDY  = ( DGDS  DGDX1  DGDX2  DGDP  DGDR  0   (DGDE3) )
!               ( DLDS  DLDX1  DLDX2  DLDP  DLDR  0   (DLDE3) )
!               ( DJDS  DJDX1  DJDX2  DJDP  DJDR  0   (DJDE3) )
!               ( DKDS  DKDX1  DKDX2  DKDP  DKDR  0   (DKDE3) )
!               ( DRDS  DRDX1  DRDX2  DRDP  DRDR  0   (DRDE3) )
!               ( 0     0      0      0     0     1   (0)     )
!               ((DQDS)(DQDX1)(DQDX2)(DQDP)(DQDR)(0)  (DQDE3) )
!                                                     ( SI IOPTIO = 0 )
!
!                   DY  = ( DSIG DX1 DX2 DP DR DQ DXXI (DEPS3))
!                   Y   = ( SIG  X1  X2  P  R  Q  XXI  (EPS3) )
!
!       DRDY  = ( DGDS  DGDX1  DGDX2  DGDP  DGDR  DGDQ  DGDXXI (DGDE3) )
!               ( DLDS  DLDX1  DLDX2  DLDP  DLDR  DLDQ  DLDXXI (DLDE3) )
!               ( DJDS  DJDX1  DJDX2  DJDP  DJDR  DJDQ  DJDXXI (DJDE3) )
!               ( DKDS  DKDX1  DKDX2  DKDP  DKDR  DKDQ  DKDXXI (DKDE3) )
!               ( DRDS  DRDX1  DRDX2  DRDP  DRDR  DRDQ  DRDXXI (DRDE3) )
!               ( DTDS  DTDX1  DTDX2  DTDP  DTDR  DTDQ  DTDXXI (DTDE3) )
!               ( DXIDS DXIDX1 DXIDX2 DXIDP DXIDR DXIDQ DXIDXI(DXIDE3))
!               ((DQDS)(DQDX1)(DQDX2)(DQDP)(DQDR)(DQDQ)(DQDXXI)(DQDE3) )
!                                                     ( SI IOPTIO = 2 )
!
!       IN  MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           YF     :  VARIABLES A T + DT = ( SIGF  X1F X2F PF RF QF(..))
!           DY     :  SOLUTION ESSAI     = ( DSIG  DX1 DX2 DP DR DQ(..))
!           NMOD   :  DIMENSION DECLAREE DRDY
!           EPSD   :  DEFORMATION A T
!           DEPS   :  INCREMENT DE DEFORMATION
!       OUT DRDY   :  JACOBIEN DU SYSTEME NON LINEAIRE
!       ----------------------------------------------------------------
!
#include "asterfort/chbfs.h"
#include "asterfort/chbfss.h"
#include "asterfort/chbfsx.h"
#include "asterfort/cvmcvx.h"
#include "asterfort/lcicma.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcprte.h"
    integer(kind=8) :: ndt, ndi, nmat, nmod
    integer(kind=8) :: ioptio, idnr, nopt
!
    real(kind=8) :: un, zero, d23, d13, mun
    parameter(un=1.d0)
    parameter(mun=-1.d0)
    parameter(zero=0.d0)
    parameter(d23=2.d0/3.d0)
    parameter(d13=-1.d0/3.d0)
!
    real(kind=8) :: hook(6, 6), ddfdds(6, 6), ddfdsx(6, 6), i6(6, 6)
    real(kind=8) :: deps(6), epsd(6), fkooh(6, 6), id(6, 6)
    real(kind=8) :: dfds(6)
    real(kind=8) :: sig(6), dsig(6)
    real(kind=8) :: x1(6), dx1(6), x2(6), dx2(6)
    real(kind=8) :: xxi(6), dxxi(6), p, dp
    real(kind=8) :: r, q
    real(kind=8) :: yf(*), dy(*), drdy(nmod, nmod)
!
    real(kind=8) :: dgds(6, 6), dgdx1(6, 6), dgdx2(6, 6), dgdr(6)
    real(kind=8) :: dlds(6, 6), dldx1(6, 6), dldx2(6, 6), dldr(6)
    real(kind=8) :: djds(6, 6), djdx1(6, 6), djdx2(6, 6), djdr(6)
    real(kind=8) :: dkds(6), dkdx1(6), dkdx2(6), dkdr
    real(kind=8) :: drds(6), drdx1(6), drdx2(6), drdr
    real(kind=8) :: dtds(6), dtdx1(6), dtdx2(6), dtdr
    real(kind=8) :: dxids(6, 6), dxidx1(6, 6), dxidx2(6, 6), dxidr(6)
    real(kind=8) :: dqds(6), dqdx1(6), dqdx2(6), dqdr
!
    real(kind=8) :: dgdq(6), dgdp(6), dgdxxi(6, 6), dgde3(6)
    real(kind=8) :: dldq(6), dldp(6), dldxxi(6, 6), dlde3(6)
    real(kind=8) :: djdq(6), djdp(6), djdxxi(6, 6), djde3(6)
    real(kind=8) :: dkdq, dkdp, dkdxxi(6), dkde3
    real(kind=8) :: drdq, drdp, drdxxi(6), drde3
    real(kind=8) :: dtdq, dtdp, dtdxxi(6), dtde3
    real(kind=8) :: dxidq(6), dxidp(6), dxidxi(6, 6), dxide3(6)
    real(kind=8) :: dqdq, dqdp, dqdxxi(6), dqde3
!
    real(kind=8) :: materf(nmat, 2), seuil
    real(kind=8) :: timed, timef, dt
!
    real(kind=8) :: k0, ak, ar, n, alp, ww
    real(kind=8) :: b, mr, gr, mu, qm, q0
    real(kind=8) :: qr0, eta, ai
    real(kind=8) :: m1, d1, gx1, g10, c1
    real(kind=8) :: m2, d2, gx2, g20, c2
    real(kind=8) :: ccin, xx, yy, zz
    real(kind=8) :: grq, qr
    real(kind=8) :: c1d, c2d, difc1, difc2
    real(kind=8) :: vtmp(6), vtmp1(6), epsp(6), epxi(6)
    real(kind=8) :: dede3(6), mtmp(6, 6), mtmp1(6, 6)
    real(kind=8) :: x1df, jx1, x2df, jx2, dcin
    real(kind=8) :: jepxi, epxino(6), nnet
    integer(kind=8) :: n1, n2, n3, n4, n5, n6, n7, n8
!
    character(len=8) :: mod
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/opti/ioptio, idnr
    common/coed/c1d, c2d
!       ----------------------------------------------------------------
    data dede3/zero, zero, mun, zero, zero, zero/
    data i6/un, zero, zero, zero, zero, zero,&
     &                   zero, un, zero, zero, zero, zero,&
     &                   zero, zero, un, zero, zero, zero,&
     &                   zero, zero, zero, un, zero, zero,&
     &                   zero, zero, zero, zero, un, zero,&
     &                   zero, zero, zero, zero, zero, un/
    data id/d23, d13, d13, zero, zero, zero,&
     &                    d13, d23, d13, zero, zero, zero,&
     &                    d13, d13, d23, zero, zero, zero,&
     &                    zero, zero, zero, un, zero, zero,&
     &                    zero, zero, zero, zero, un, zero,&
     &                    zero, zero, zero, zero, zero, un/
!
    sig(1:ndt) = yf(1:ndt)
    x1(1:ndt) = yf(ndt+1:ndt+ndt)
    x2(1:ndt) = yf(2*ndt+1:2*ndt+ndt)
    p = yf(3*ndt+1)
    r = yf(3*ndt+2)
    q = yf(3*ndt+3)
    dsig(1:ndt) = dy(1:ndt)
    dx1(1:ndt) = dy(ndt+1:ndt+ndt)
    dx2(1:ndt) = dy(2*ndt+1:2*ndt+ndt)
    dp = dy(3*ndt+1)
!
    k0 = materf(1, 2)
    ak = materf(2, 2)
    ar = materf(3, 2)
    n = materf(5, 2)
    alp = materf(6, 2)
    b = materf(7, 2)
    mr = materf(8, 2)
    gr = materf(9, 2)
    mu = materf(10, 2)
    qm = materf(11, 2)
    q0 = materf(12, 2)
    qr0 = materf(13, 2)
    eta = materf(14, 2)
    c1 = materf(15, 2)
    m1 = materf(16, 2)
    d1 = materf(17, 2)
    gx1 = materf(18, 2)
    g10 = materf(19, 2)
    c2 = materf(20, 2)
    m2 = materf(21, 2)
    d2 = materf(22, 2)
    gx2 = materf(23, 2)
    g20 = materf(24, 2)
    ai = materf(25, 2)
!
    nopt = 0
    if (ioptio .eq. 2) nopt = idnr
!
    call lcopli('ISOTROPE', mod, materf(1, 1), hook)
    call chbfs(sig, x1, x2, dfds)
    call chbfss(sig, x1, x2, id, ddfdds)
    call chbfsx(sig, x1, x2, i6, ddfdsx)
!
    call cvmcvx(nmat, materf, sig, yf(ndt+1), seuil)
    if (seuil .lt. 0.d0) seuil = 0.d0
    ccin = ai+(1.d0-ai)*exp(-b*p)
    dcin = b*(ai-1.d0)*exp(-b*p)
    dt = timef-timed
!
!       ----------------------------------------------------------------
!       CALCUL DU JACOBIEN DU SYSTEME ( SIG  X1  X2  P  R  (EPS3) )
!       ----------------------------------------------------------------
!
! - DGDS(T+DT)
    dgds(1:ndt, 1:ndt) = dp*matmul(hook(1:ndt, 1:ndt), ddfdds(1:ndt, 1:ndt))
    dgds(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)+dgds(1:ndt, 1:ndt)
!
! - DGDX1(T+DT)
    dgdx1(1:ndt, 1:ndt) = dp*matmul(hook(1:ndt, 1:ndt), ddfdsx(1:ndt, 1:ndt))
!
! - DGDX2(T+DT)
    dgdx2(1:ndt, 1:ndt) = dgdx1(1:ndt, 1:ndt)
!
! - DGDP(T+DT)
    if (seuil .lt. 0.d0) then
        dgdp(:) = 0.d0
    else
        dgdp(1:ndt) = matmul(hook(1:ndt, 1:ndt), dfds(1:ndt))
    end if
!
! - DGDR(T+DT)
    dgdr(:) = 0.d0
!
! - DGDQ(T+DT)
    dgdq(:) = 0.d0
!
!
! - DLDS(T+DT)
    jx1 = dot_product(x1(1:ndt), x1(1:ndt))
    jx1 = sqrt(jx1*3.d0/2.d0)
    xx = c1*dp*2.d0/3.d0
    yy = g10*ccin*(1.d0-d1)*dp*2.d0/3.d0
    if (jx1 .le. 0.d0) then
        zz = 1.d0+dp*g10*ccin*d1
        ww = 0.d0
    else
        zz = 1.d0+dp*g10*ccin*d1+dt*gx1*jx1**(m1-1.d0)
        ww = gx1*dt*(m1-1.d0)*jx1**(m1-3.d0)*3.d0/2.d0
    end if
    vtmp(1:ndt) = matmul(ddfdds(1:ndt, 1:ndt), x1(1:ndt))
    call lcprte(vtmp, dfds, mtmp)
    x1df = dot_product(x1(1:ndt), dfds(1:ndt))
    mtmp1(1:ndt, 1:ndt) = x1df*ddfdds(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    dlds(1:ndt, 1:ndt) = yy*mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = xx*ddfdds(1:ndt, 1:ndt)
    dlds(1:ndt, 1:ndt) = dlds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
!
! - DLDX1(T+DT)
    vtmp(1:ndt) = matmul(ddfdsx(1:ndt, 1:ndt), x1(1:ndt))
    call lcprte(vtmp, dfds, mtmp1)
    mtmp(1:ndt, 1:ndt) = x1df*ddfdsx(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    call lcprte(dfds, dfds, mtmp)
    mtmp(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    dldx1(1:ndt, 1:ndt) = yy*mtmp(1:ndt, 1:ndt)
    call lcprte(x1, x1, mtmp)
    mtmp(1:ndt, 1:ndt) = ww*mtmp(1:ndt, 1:ndt)
    dldx1(1:ndt, 1:ndt) = dldx1(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = zz*i6(1:ndt, 1:ndt)
    dldx1(1:ndt, 1:ndt) = dldx1(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = xx*ddfdsx(1:ndt, 1:ndt)
    dldx1(1:ndt, 1:ndt) = dldx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
!
! - DLDX2(T+DT)
    dldx2(1:ndt, 1:ndt) = yy*mtmp1(1:ndt, 1:ndt)
    dldx2(1:ndt, 1:ndt) = dldx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
!
! -- CAS ANISOTHERME
    if (c1 .ne. 0.d0) then
        difc1 = (c1-c1d)/c1
        mtmp(1:ndt, 1:ndt) = difc1*i6(1:ndt, 1:ndt)
        dldx1(1:ndt, 1:ndt) = dldx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    end if
    if (c2 .ne. 0.d0) then
        difc2 = (c2-c2d)/c2
        mtmp(1:ndt, 1:ndt) = difc2*i6(1:ndt, 1:ndt)
        dldx2(1:ndt, 1:ndt) = dldx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    end if
!
! - DLDP(T+DT)
    yy = g10*(ccin+dcin*dp)*d1
    zz = g10*(ccin+dcin*dp)*(1.d0-d1)*2.d0/3.d0
    xx = x1df*zz-c1*2.d0/3.d0
    vtmp(1:ndt) = xx*dfds(1:ndt)
    dldp(1:ndt) = yy*x1(1:ndt)
    dldp(1:ndt) = dldp(1:ndt)+vtmp(1:ndt)
!
! - DLDR(T+DT)
    dldr(:) = 0.d0
!
! - DLDQ(T+DT)
    dldq(:) = 0.d0
!
!
! - DJDS(T+DT)
    jx2 = dot_product(x2(1:ndt), x2(1:ndt))
    jx2 = sqrt(jx2*3.d0/2.d0)
    xx = c2*dp*2.d0/3.d0
    yy = g20*ccin*(1.d0-d2)*dp*2.d0/3.d0
    if (jx2 .le. 0.d0) then
        zz = 1.d0+dp*g20*ccin*d2
        ww = 0.d0
    else
        zz = 1.d0+dp*g20*ccin*d2+dt*gx2*jx2**(m2-1.d0)
        ww = gx2*dt*(m2-1.d0)*jx2**(m2-3.d0)*3.d0/2.d0
    end if
    vtmp(1:ndt) = matmul(ddfdds(1:ndt, 1:ndt), x2(1:ndt))
    call lcprte(vtmp, dfds, mtmp)
    x2df = dot_product(x2(1:ndt), dfds(1:ndt))
    mtmp1(1:ndt, 1:ndt) = x2df*ddfdds(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    djds(1:ndt, 1:ndt) = yy*mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = xx*ddfdds(1:ndt, 1:ndt)
    djds(1:ndt, 1:ndt) = djds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
!
! - DJDX2(T+DT)
    vtmp(1:ndt) = matmul(ddfdsx(1:ndt, 1:ndt), x2(1:ndt))
    call lcprte(vtmp, dfds, mtmp1)
    mtmp(1:ndt, 1:ndt) = x2df*ddfdsx(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    call lcprte(dfds, dfds, mtmp)
    mtmp(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    djdx2(1:ndt, 1:ndt) = yy*mtmp(1:ndt, 1:ndt)
    call lcprte(x2, x2, mtmp)
    mtmp(1:ndt, 1:ndt) = ww*mtmp(1:ndt, 1:ndt)
    djdx2(1:ndt, 1:ndt) = djdx2(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = zz*i6(1:ndt, 1:ndt)
    djdx2(1:ndt, 1:ndt) = djdx2(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = xx*ddfdsx(1:ndt, 1:ndt)
    djdx2(1:ndt, 1:ndt) = djdx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
!
! - DJDX1(T+DT)
    djdx1(1:ndt, 1:ndt) = yy*mtmp1(1:ndt, 1:ndt)
    djdx1(1:ndt, 1:ndt) = djdx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
!
! - DJDP(T+DT)
    yy = g20*(ccin+dcin*dp)*d2
    zz = g20*(ccin+dcin*dp)*(1.d0-d2)*2.d0/3.d0
    xx = x2df*zz-c2*2.d0/3.d0
    vtmp(1:ndt) = xx*dfds(1:ndt)
    djdp(1:ndt) = yy*x2(1:ndt)
    djdp(1:ndt) = djdp(1:ndt)+vtmp(1:ndt)
!
! - DJDR(T+DT)
    djdr(:) = 0.d0
!
! - DJDQ(T+DT)
    djdq(:) = 0.d0
!
!
! - DKDS(T+DT)
    xx = seuil/(k0+ak*r)
    if (xx .lt. 0.d0) xx = 0.d0
    zz = dt*( &
         (xx**(n-1.d0))*(n+alp*(n+1)*xx**(n+1)))*exp(alp*(xx**(n+1)))/(k0+ak*r)
    dkds(1:ndt) = (-zz)*dfds(1:ndt)
!
! - DKDX1(T+DT)
    dkdx1(1:ndt) = zz*dfds(1:ndt)
!
! - DKDX2(T+DT)
    dkdx2(1:ndt) = dkdx1(1:ndt)
!
! - DKDP(T+DT)
    dkdp = 1.d0
!
! - DKDR(T+DT)
!       DKDR = ZZ * (AR * K0 + SEUIL - AK * K) / (K0 + AK * R)
    dkdr = zz*(ar*(k0+ak*r)+ak*seuil)/(k0+ak*r)
!
! - DKDQ(T+DT)
    dkdq = 0.d0
!
!
! - DRDS(T+DT)
    drds(:) = 0.d0
!
! - DRDX1(T+DT)
    drdx1(:) = 0.d0
!
! - DRDX2(T+DT)
    drdx2(:) = 0.d0
!
! - DRDP(T+DT)
    grq = q0+(qm-q0)*(1.d0-exp(-2.d0*mu*q))
    qr = grq-qr0*(1.d0-((qm-grq)/qm)**2)
    drdp = b*(r-grq)
!
! - DRDR(T+DT)
    drdr = 1.d0+b*dp+gr*dt*mr*(abs(qr-r))**(mr-1.d0)
!
! - DRDQ(T+DT)
    drdq = 0.d0
!
!
! - DTDS(T+DT)
    dtds(:) = 0.d0
!
! - DTDX1(T+DT)
    dtdx1(:) = 0.d0
!
! - DTDX2(T+DT)
    dtdx2(:) = 0.d0
!
! - DTDP(T+DT)
    dtdp = 0.d0
!
! - DTDR(T+DT)
    dtdr = 0.d0
!
! - DTDQ(T+DT)
    dtdq = 1.d0
!
! - CONTRAINTES PLANES -------------------------------------------------
!
    if (mod(1:6) .eq. 'C_PLAN') then
!
! - DGDE3(T+DT)
        dgde3(1:ndt) = matmul(hook(1:ndt, 1:ndt), dede3(1:ndt))
!
! - DLDE3(T+DT)
        dlde3(:) = 0.d0
!
! - DJDE3(T+DT)
        djde3(:) = 0.d0
!
! - DKDE3(T+DT)
        dkde3 = 0.d0
!
! - DRDE3(T+DT)
        drde3 = 0.d0
!
! - DTDE3(T+DT)
        dtde3 = 0.d0
!
! - DQDE3(T+DT)
        dqde3 = hook(3, 3)
!
! - DQDS (T+DT)
        dqds(1) = -dp*( &
                  hook(3, 3)*ddfdds(3, 1)+hook(3, 1)*ddfdds(1, 1)+hook(3, 2)*ddfdds(2, 1)+hook(3, 4&
                  &)*ddfdds(4, 1) &
                  )
        dqds(2) = -dp*( &
                  hook(3, 3)*ddfdds(3, 2)+hook(3, 1)*ddfdds(1, 2)+hook(3, 2)*ddfdds(2, 2)+hook(3, 4&
                  &)*ddfdds(4, 2) &
                  )
        dqds(3) = -dp*( &
                  hook(3, 3)*ddfdds(3, 3)+hook(3, 1)*ddfdds(1, 3)+hook(3, 2)*ddfdds(2, 3)+hook(3, 4&
                  &)*ddfdds(4, 3) &
                  )
        dqds(4) = -dp*( &
                  hook(3, 3)*ddfdds(3, 4)+hook(3, 1)*ddfdds(1, 4)+hook(3, 2)*ddfdds(2, 4)+hook(3, 4&
                  &)*ddfdds(4, 4) &
                  )
!
! - DQDX1 (T+DT)
        dqdx1(1) = -dp*(hook(3, 3)*ddfdsx(3, 1)+hook(3, 1)*ddfdsx(1, 1)+ &
                        hook(3, 2)*ddfdsx(2, 1)+hook(3, 4)*ddfdsx(4, 1))
        dqdx1(2) = -dp*(hook(3, 3)*ddfdsx(3, 2)+hook(3, 1)*ddfdsx(1, 2)+ &
                        hook(3, 2)*ddfdsx(2, 2)+hook(3, 4)*ddfdsx(4, 2))
        dqdx1(3) = -dp*(hook(3, 3)*ddfdsx(3, 3)+hook(3, 1)*ddfdsx(1, 3)+ &
                        hook(3, 2)*ddfdsx(2, 3)+hook(3, 4)*ddfdsx(4, 3))
        dqdx1(4) = -dp*(hook(3, 3)*ddfdsx(3, 4)+hook(3, 1)*ddfdsx(1, 4)+ &
                        hook(3, 2)*ddfdsx(2, 4)+hook(3, 4)*ddfdsx(4, 4))
!
! - DQDX2 (T+DT)
        dqdx2(1:ndt) = dqdx1(1:ndt)
!
! - DQDP (T+DT)
        dqdp = -hook(3, 1)*dfds(1)-hook(3, 2)*dfds(2)-hook(3, 3)*dfds(3)-hook(3, 4)*dfds(4)
!
! - DQDR (T+DT)
        dqdr = 0.d0
!
! - DQDQ (T+DT)
        dqdq = 0.d0
!
    end if
!
!
! - ASSEMBLAGE ---------------------------------------------------------
!
! - DRDY (T+DT) = ( DGDS  DGDX1  DGDX2  DGDP  DGDR  0   (DGDE3) )
!                 ( DLDS  DLDX1  DLDX2  DLDP  DLDR  0   (DLDE3) )
!                 ( DJDS  DJDX1  DJDX2  DJDP  DJDR  0   (DJDE3) )
!                 ( DKDS  DKDX1  DKDX2  DKDP  DKDR  0   (DKDE3) )
!                 ( DRDS  DRDX1  DRDX2  DRDP  DRDR  0   (DRDE3) )
!                 ( 0     0      0      0     0     1   (0)     )
!                 ((DQDS)(DQDX1)(DQDX2)(DQDP)(DQDR)(0)  (DQDE3) )
!
!
    n1 = 1
    n2 = ndt+1
    n3 = 2*ndt+1
    n4 = 3*ndt+1
    n5 = 3*ndt+2
    n6 = 3*ndt+3
    n8 = 3*ndt+4+nopt
!
    call lcicma(dgds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n1, n1)
    call lcicma(dgdx1, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n1, n2)
    call lcicma(dgdx2, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n1, n3)
    call lcicma(dgdp, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n1, n4)
    call lcicma(dgdr, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n1, n5)
    call lcicma(dgdq, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n1, n6)
!
    call lcicma(dlds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n2, n1)
    call lcicma(dldx1, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n2, n2)
    call lcicma(dldx2, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n2, n3)
    call lcicma(dldp, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n2, n4)
    call lcicma(dldr, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n2, n5)
    call lcicma(dldq, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n2, n6)
!
    call lcicma(djds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n3, n1)
    call lcicma(djdx1, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n3, n2)
    call lcicma(djdx2, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                n3, n3)
    call lcicma(djdp, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n3, n4)
    call lcicma(djdr, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n3, n5)
    call lcicma(djdq, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                n3, n6)
!
    call lcicma(dkds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n4, n1)
    call lcicma(dkdx1, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n4, n2)
    call lcicma(dkdx2, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n4, n3)
    drdy(n4, n4) = dkdp
    drdy(n4, n5) = dkdr
    drdy(n4, n6) = dkdq
!
    call lcicma(drds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n5, n1)
    call lcicma(drdx1, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n5, n2)
    call lcicma(drdx2, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n5, n3)
    drdy(n5, n4) = drdp
    drdy(n5, n5) = drdr
    drdy(n5, n6) = drdq
!
    call lcicma(dtds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n6, n1)
    call lcicma(dtdx1, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n6, n2)
    call lcicma(dtdx2, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                n6, n3)
    drdy(n6, n4) = dtdp
    drdy(n6, n5) = dtdr
    drdy(n6, n6) = dtdq
!
    if (mod(1:6) .eq. 'C_PLAN') then
!
        call lcicma(dgde3, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n1, n8)
        call lcicma(dlde3, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n2, n8)
        call lcicma(djde3, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n3, n8)
        drdy(n4, n8) = dkde3
        drdy(n5, n8) = drde3
        drdy(n6, n8) = dtde3
!
        call lcicma(dqds, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n8, n1)
        call lcicma(dqdx1, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n8, n2)
        call lcicma(dqdx2, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n8, n3)
        drdy(n8, n4) = dqdp
        drdy(n8, n5) = dqdr
        drdy(n8, n6) = dqdq
        drdy(n8, n8) = dqde3
    end if
!
!       ----------------------------------------------------------------
!       CALCUL DU JACOBIEN DU SYSTEME (SIG  X1  X2  P  R  Q XXI (EPS3))
!       ----------------------------------------------------------------
!
    if (ioptio .eq. 2) then
!
        xxi(1:ndt) = yf(3*ndt+4:3*ndt+4-1+ndt)
        dxxi(1:ndt) = dy(3*ndt+4:3*ndt+4-1+ndt)
!
! - DGDQ(T+DT)
        dgdq(:) = 0.d0
!
! - DLDQ(T+DT)
        dldq(:) = 0.d0
!
! - DJDQ(T+DT)
        djdq(:) = 0.d0
!
! - DKDQ(T+DT)
        dkdq = 0.d0
!
! - DRDQ(T+DT)
        xx = gr*mr*((abs(qr-r))**(mr-1.d0))*(1.d0-2.d0*(qm-grq)*qr0/(qm**2))
        drdq = (q0-qm)*2.d0*mu*exp(-2.d0*mu*q)*(b*dp+xx*dt)
!
! - DTDQ(T+DT)
        dtdq = 1.d0
!
!
! - DGDXXI(T+DT)
        dgdxxi(:, :) = 0.d0
!
! - DLDXXI(T+DT)
        dldxxi(:, :) = 0.d0
!
! - DJDXXI(T+DT)
        djdxxi(:, :) = 0.d0
!
! - DKDXXI(T+DT)
        dkdxxi(:) = 0.d0
!
! - DRDXXI(T+DT)
        drdxxi(:) = 0.d0
!
!
!
! ---- EPSP
        call lcopil('ISOTROPE', mod, materf(1, 1), fkooh)
        vtmp(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sig(1:ndt))
        epsp(1:ndt) = epsd(1:ndt)+deps(1:ndt)
        epsp(1:ndt) = epsp(1:ndt)-vtmp(1:ndt)
! ---- JEPXI
        epxi(1:ndt) = epsp(1:ndt)-xxi(1:ndt)
        xx = dot_product(epxi(1:ndt), epxi(1:ndt))
        jepxi = sqrt(xx*3.d0/2.d0)
!
! --- H(F)=SEUIL2
!
!           SEUIL2 = 2.D0/3.D0 * JEPXI - Q
! --- NNET
!
        if (jepxi .eq. 0.d0) then
            nnet = 0.d0
        else
            epxino(1:ndt) = (1.d0/jepxi)*epxi(1:ndt)
            nnet = dot_product(dfds(1:ndt), epxino(1:ndt))
        end if
!
! -MEMORISATION
!
        if (jepxi .gt. 0.d0) then
!
            zz = -eta/jepxi
            xx = zz*dp
            vtmp(1:ndt) = epsp(1:ndt)-xxi(1:ndt)
            vtmp1(1:ndt) = dp*dfds(1:ndt)
            vtmp(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
            yy = -dp*nnet*(3.d0/2.d0)/jepxi
            vtmp1(1:ndt) = yy*epxi(1:ndt)
            vtmp(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
            vtmp1(1:ndt) = xx*vtmp(1:ndt)
!
! - DTDS(T+DT)
!
            dtds(1:ndt) = matmul(ddfdds(1:ndt, 1:ndt), vtmp1(1:ndt))
!
! - DTDX1(T+DT)
            dtdx1(1:ndt) = matmul(ddfdsx(1:ndt, 1:ndt), vtmp1(1:ndt))
!
! - DTDX2(T+DT)
            dtdx2(1:ndt) = dtdx1(1:ndt)
!
! - DTDP(T+DT)
!
            vtmp(1:ndt) = zz*vtmp(1:ndt)
            dtdp = dot_product(dfds(1:ndt), vtmp(1:ndt))
!
! - DTDXXI(T+DT)
            yy = -nnet*(3.d0/2.d0)/jepxi
            vtmp(1:ndt) = yy*epxi(1:ndt)
            vtmp(1:ndt) = vtmp(1:ndt)+dfds(1:ndt)
            dtdxxi(1:ndt) = (-xx)*vtmp(1:ndt)
!
!
! - DXIDS(T+DT)
            zz = 3.d0/2.d0*(1.d0-eta)*(1.d0-dp*nnet/jepxi*3.d0)
            xx = 3.d0/2.d0*(1.d0-eta)*dp/jepxi
            vtmp(1:ndt) = (-zz*dp)*epxi(1:ndt)
            vtmp1(1:ndt) = (-xx*dp)*dfds(1:ndt)
            vtmp(1:ndt) = vtmp(1:ndt)+vtmp1(1:ndt)
            vtmp1(1:ndt) = matmul(ddfdds(1:ndt, 1:ndt), vtmp(1:ndt))
            call lcprte(vtmp1, epxi, mtmp)
            mtmp1(1:ndt, 1:ndt) = (-xx*dp*nnet)*ddfdds(1:ndt, 1:ndt)
            dxids(1:ndt, 1:ndt) = mtmp1(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
!
! - DXIDX1(T+DT)
            vtmp1(1:ndt) = matmul(ddfdsx(1:ndt, 1:ndt), vtmp(1:ndt))
            call lcprte(vtmp1, epxi, mtmp)
            mtmp1(1:ndt, 1:ndt) = (-xx*dp*nnet)*ddfdsx(1:ndt, 1:ndt)
            dxidx1(1:ndt, 1:ndt) = mtmp1(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
!
! - DXIDX2(T+DT)
            dxidx2(1:ndt, 1:ndt) = dxidx1(1:ndt, 1:ndt)
!
! - DXIDP(T+DT)
            yy = dot_product(dfds(1:ndt), dfds(1:ndt))
            vtmp(1:ndt) = (zz*nnet+xx*yy)*epxi(1:ndt)
            vtmp1(1:ndt) = (-xx*nnet)*dfds(1:ndt)
            dxidp(1:ndt) = vtmp1(1:ndt)-vtmp(1:ndt)
!
! - DXIDXI(T+DT)
            mtmp(1:ndt, 1:ndt) = (1.d0+xx*nnet)*i6(1:ndt, 1:ndt)
            vtmp(1:ndt) = (3.d0*xx*nnet)*epxi(1:ndt)
            vtmp1(1:ndt) = xx*dfds(1:ndt)
            vtmp(1:ndt) = vtmp1(1:ndt)-vtmp(1:ndt)
            call lcprte(vtmp, epxi, mtmp1)
            dxidxi(1:ndt, 1:ndt) = mtmp1(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
!
        else
            dgdxxi(:, :) = 0.d0
!
! - DTDS(T+DT)
            dtds(:) = 0.d0
!
! - DTDX1(T+DT)
            dtdx1(:) = 0.d0
!
! - DTDX2(T+DT)
            dtdx2(:) = 0.d0
!
! - DTDP(T+DT)
            dtdp = 0.d0
!
! - DTDXXI(T+DT)
            dtdxxi(:) = 0.d0
!
! - DXIDS(T+DT)
            dxids(:, :) = 0.d0
!
! - DXIDX1(T+DT)
            dxidx1(:, :) = 0.d0
!
! - DXIDX2(T+DT)
            dxidx2(:, :) = 0.d0
!
! - DXIDP(T+DT)
            dxidp(:) = 0.d0
!
! - DXIDXI(T+DT)
            dxidxi(:, :) = 0.d0
!
        end if
!
! - DTDR(T+DT)
        dtdr = 0.d0
!
! - DTDQ(T+DT)
        dtdq = 1.d0
!
! - DXIDR(T+DT)
        dxidr(:) = 0.d0
!
! - DXIDQ(T+DT)
        dxidq(:) = 0.d0
!
!
        if (mod(1:6) .eq. 'C_PLAN') then
!
! - DQDXXI(T+DT)
            dqdxxi(:) = 0.d0
!
! - DXIDE3(T+DT)
            dxide3(:) = 0.d0
!
! - DQDQ(T+DT)
            dqdq = 0.d0
!
! - DTDE3(T+DT)
            dtde3 = 0.d0
        end if
!
!
! - ASSEMBLAGE ---------------------------------------------------------
!
!       DRDY  = (                                 DGDQ  DGDXXI (DGDE3) )
!               (         (...)                   DLDQ  DLDXXI (DLDE3) )
!               (         (...)                   DJDQ  DJDXXI (DJDE3) )
!               (                                 DKDQ  DKDXXI (DKDE3) )
!               (                                 DRDQ  DRDXXI (DRDE3) )
!               ( DTDS  DTDX1  DTDX2  DTDP  DTDR  DTDQ  DTDXXI (DTDE3) )
!               ( DXIDS DXIDX1 DXIDX2 DXIDP DXIDR DXIDQ DXIDXI(DXIDE3))
!               ((DQDS)(DQDX1)(DQDX2)(DQDP)(DQDR)(DQDQ)(DQDXXI)(DQDE3) )
!
!
        n7 = 3*ndt+4
!
        call lcicma(dgdq, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n1, n6)
        call lcicma(dgdxxi, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n1, n7)
!
        call lcicma(dldq, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n2, n6)
        call lcicma(dldxxi, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n2, n7)
!
        call lcicma(djdq, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n3, n6)
        call lcicma(djdxxi, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n3, n7)
!
        drdy(n4, n6) = dkdq
        call lcicma(dkdxxi, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n4, n7)
!
        drdy(n5, n6) = drdq
        call lcicma(drdxxi, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n5, n7)
!
        call lcicma(dtds, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n6, n1)
        call lcicma(dtdx1, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n6, n2)
        call lcicma(dtdx2, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n6, n3)
        drdy(n6, n4) = dtdp
        drdy(n6, n5) = dtdr
        drdy(n6, n6) = dtdq
        call lcicma(dtdxxi, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n6, n7)
!
        call lcicma(dxids, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n1)
        call lcicma(dxidx1, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n2)
        call lcicma(dxidx2, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n3)
        call lcicma(dxidp, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n4)
        call lcicma(dxidr, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n5)
        call lcicma(dxidq, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n6)
        call lcicma(dxidxi, 6, 6, ndt, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    n7, n7)
!
        if (mod(1:6) .eq. 'C_PLAN') then
            drdy(n8, n6) = dqdq
            call lcicma(dqdxxi, 1, 6, 1, ndt, &
                        1, 1, drdy, nmod, nmod, &
                        n8, n7)
            call lcicma(dxide3, 6, 1, ndt, 1, &
                        1, 1, drdy, nmod, nmod, &
                        n7, n8)
            drdy(n6, n8) = dtde3
        end if
!
    end if
!
end subroutine
