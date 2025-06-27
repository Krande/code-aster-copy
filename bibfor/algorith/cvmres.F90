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

subroutine cvmres(mod, nmat, materd, materf, timed, &
                  timef, yd, yf, epsd, deps, &
                  dy, res)
    implicit none
!     VISCOCHABOCHE :
!            CALCUL DES TERMES DU SYSTEME NL A RESOUDRE = RES(DY)
!                   DY  = ( DSIG DX1 DX2 DP DR 0  (DEPS3))
!                   Y   = ( SIG  X1  X2  P  R  Q  (EPS3) )
!                   RES = ( GF   LF  JF  KF RF 0  (F)    )
!                                                      (SI IOPTIO = 0)
!
!                   DY =  ( DSIG DX1 DX2 DP DR DQ DXXI (DEPS3))
!                   Y  =  ( SIG  X1  X2  P  R  Q  XXI  (EPS3) )
!                   RES = ( GF   LF  JF  KF RF TF XIF  (FF)   )
!                                                      (SI IOPTIO = 2)
!
!       IN  MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERD :  COEFFICIENTS MATERIAU A T
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           YD     :  VARIABLES A T      = ( SIGD  X1D X2D PD RD QD(..))
!           YF     :  VARIABLES A T + DT = ( SIGF  X1F X2F PF RF QF(..))
!           DY     :  SOLUTION ESSAI     = ( DSIG  DX1 DX2 DP DR DQ(..))
!           EPSD   :  DEFORMATION A T
!           DEPS   :  INCREMENT DE DEFORMATION
!       OUT RES    :  SYSTEME NL A T + DT
!       ----------------------------------------------------------------
#include "asterfort/chbfs.h"
#include "asterfort/cvmcvx.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
    integer(kind=8) :: ndt, ndi, nmat
    integer(kind=8) :: ioptio, idnr, nopt
!
    real(kind=8) :: hookf(6, 6), dkooh(6, 6), fkooh(6, 6)
    real(kind=8) :: sigf(6), dsig(6), sigd(6), dfds(6)
    real(kind=8) :: deps(6), depsp(6), depse(6)
    real(kind=8) :: epsed(6), epsef(6), epsd(6)
    real(kind=8) :: x1(6), dx1(6), x2(6), dx2(6)
    real(kind=8) :: xxi(6), dxxi(6), p, dp
    real(kind=8) :: r, dr, q, dq
    real(kind=8) :: gf(6), lf(6), jf(6), kf
    real(kind=8) :: rf, tf, xif(6), ff
    real(kind=8) :: res(*), dy(*), yd(*), yf(*)
    real(kind=8) :: timed, timef, dt
!
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: seuil
!
    real(kind=8) :: k0, ak, n, alp
    real(kind=8) :: b, mr, gr, mu, qm, q0
    real(kind=8) :: qr0, eta, ai
    real(kind=8) :: m1, d1, gx1, g10, c1, c1d
    real(kind=8) :: m2, d2, gx2, g20, c2, c2d
    real(kind=8) :: difc1, difc2
    real(kind=8) :: ccin, xx, yy, zz, sgn
    real(kind=8) :: grq, qr
    real(kind=8) :: vtmp(6), vtmp1(6), epsp(6), epxi(6)
!
    character(len=8) :: mod
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/opti/ioptio, idnr
!       ----------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    sigd(1:ndt) = yd(1:ndt)
    sigf(1:ndt) = yf(1:ndt)
    x1(1:ndt) = yf(ndt+1:ndt+ndt)
    x2(1:ndt) = yf(2*ndt+1:2*ndt+ndt)
    p = yf(3*ndt+1)
    r = yf(3*ndt+2)
    q = yf(3*ndt+3)
    dsig(1:ndt) = dy(1:ndt)
    dx1(1:ndt) = dy(ndt+1:ndt+ndt)
    dx2(1:ndt) = dy(2*ndt+1:2*ndt+ndt)
    dp = dy(3*ndt+1)
    dr = dy(3*ndt+2)
    dq = dy(3*ndt+3)
!
!
    k0 = materf(1, 2)
    ak = materf(2, 2)
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
    c1d = materd(15, 2)
    c2d = materd(20, 2)
!
    nopt = 0
    if (ioptio .eq. 2) nopt = idnr
!
!                   -1                     -1
! -     HOOKF, HOOKD , DFDS , EPSEF = HOOKD  SIGD + DEPS - DEPSP
!
    call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
    call lcopil('ISOTROPE', mod, materd(1, 1), dkooh)
!
    call chbfs(sigf, x1, x2, dfds)
    call cvmcvx(nmat, materf, sigf, yf(ndt+1), seuil)
    ccin = ai+(1.d0-ai)*exp(-b*p)
    dt = timef-timed
!
! -     CALCUL DU RESIDU POUR ( SIG  X1  X2  P  R  0  (EPS3))
!
! - GF  (T+DT)
!
    depsp(1:ndt) = dp*dfds(1:ndt)
    epsed(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigd(1:ndt))
    depse(1:ndt) = deps(1:ndt)-depsp(1:ndt)
    epsef(1:ndt) = epsed(1:ndt)+depse(1:ndt)
    gf(1:ndt) = matmul(hookf(1:ndt, 1:ndt), epsef(1:ndt))
    gf(1:ndt) = gf(1:ndt)-sigf(1:ndt)
!
! - LF (T+DT)
!
    zz = dot_product(x1(1:ndt), dfds(1:ndt))
    yy = dot_product(x1(1:ndt), x1(1:ndt))
    zz = zz*(1.d0-d1)*g10*ccin*dp*2.d0/3.d0
    xx = c1*dp*2.d0/3.d0-zz
    yy = gx1*dt*(sqrt(yy*3.d0/2.d0))**(m1-1.d0)+g10*ccin*d1*dp
    lf(1:ndt) = xx*dfds(1:ndt)
    vtmp(1:ndt) = yy*x1(1:ndt)
    lf(1:ndt) = lf(1:ndt)-vtmp(1:ndt)
    lf(1:ndt) = lf(1:ndt)-dx1(1:ndt)
!
    if (c1d .ne. 0.d0) then
        difc1 = (c1-c1d)/c1
        vtmp(1:ndt) = difc1*x1(1:ndt)
        lf(1:ndt) = lf(1:ndt)+vtmp(1:ndt)
    end if
!
! - JF (T+DT)
!
    zz = dot_product(x2(1:ndt), dfds(1:ndt))
    yy = dot_product(x2(1:ndt), x2(1:ndt))
    zz = zz*(1.d0-d2)*g20*ccin*dp*2.d0/3.d0
    xx = c2*dp*2.d0/3.d0-zz
    yy = gx2*dt*(sqrt(yy*3.d0/2.d0))**(m2-1.d0)+g20*ccin*d2*dp
    jf(1:ndt) = xx*dfds(1:ndt)
    vtmp(1:ndt) = yy*x2(1:ndt)
    jf(1:ndt) = jf(1:ndt)-vtmp(1:ndt)
    jf(1:ndt) = jf(1:ndt)-dx2(1:ndt)
!
! - CAS ANISOTHERME
!
    if (c2 .ne. 0.d0) then
        difc2 = (c2-c2d)/c2
        vtmp(1:ndt) = difc2*x2(1:ndt)
        jf(1:ndt) = jf(1:ndt)+vtmp(1:ndt)
    end if
!
! - KF (T+DT)
!
    zz = seuil/(k0+ak*r)
    if (zz .lt. 0.d0) zz = 0.d0
    kf = dt*(zz**n)*exp(alp*(zz**(n+1)))-dp
!
! - RF (T+DT)
!
    grq = q0+(qm-q0)*(1.d0-exp(-2.d0*mu*q))
    qr = grq-qr0*(1.d0-((qm-grq)/qm)**2)
    sgn = 1.d0
    if ((qr-r) .lt. 0.d0) sgn = -1.d0
    rf = b*(grq-r)*dp+sgn*gr*dt*(abs(qr-r))**mr-dr
!
! - FF (T+DT)
!
    if (mod(1:6) .eq. 'C_PLAN') then
        ff = -hookf(3, 3)*epsef(3)-hookf(3, 1)*epsef(1)-hookf(3, 2)*epsef(2)-hookf(3, 4&
             &)*epsef(4)
    end if
!
! - RES (T+DT) = ( GF LF JF KF RF 0 (FF) )
!
    res(1:ndt) = gf(1:ndt)
    res(ndt+1:ndt+ndt) = lf(1:ndt)
    res(2*ndt+1:2*ndt+ndt) = jf(1:ndt)
    res(3*ndt+1) = kf
    res(3*ndt+2) = rf
    res(3*ndt+3) = 0.d0
    if (mod(1:6) .eq. 'C_PLAN') res(3*ndt+4+nopt) = ff
!
!
! -     CALCUL DU RESIDU POUR ( SIG  X1  X2  P  R  Q XXI (EPS3))
    if (ioptio .eq. 2) then
        xxi(1:ndt) = yf(3*ndt+4:3*ndt+4-1+ndt)
        dxxi(1:ndt) = dy(3*ndt+4:3*ndt+4-1+ndt)
!
! - EPSP
!
        call lcopil('ISOTROPE', mod, materf(1, 1), fkooh)
        vtmp1(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigf(1:ndt))
        epsp(1:ndt) = epsd(1:ndt)+deps(1:ndt)
        epsp(1:ndt) = epsp(1:ndt)-vtmp1(1:ndt)
!
! N-ETOILE
!
        vtmp(1:ndt) = epsp(1:ndt)-xxi(1:ndt)
        xx = dot_product(vtmp(1:ndt), vtmp(1:ndt))
        xx = sqrt(xx*3.d0/2.d0)
!
! H(F)
        zz = 2.d0/3.d0*xx-q
!
        if (zz .lt. 0.d0) then
            tf = 0.d0
            xif(:) = 0.d0
        else
!
            if (xx .eq. 0.d0) then
                epxi(:) = 0.d0
            else
                epxi(1:ndt) = (1.d0/xx)*vtmp(1:ndt)
            end if
!
! N *  N-ETOILE
!
            zz = dot_product(dfds(1:ndt), epxi(1:ndt))
!
            if (zz .le. 0.d0) then
                tf = 0.d0
                xif(:) = 0.d0
            else
!
! - TF (T+DT)
!
                xx = dot_product(dfds(1:ndt), epxi(1:ndt))
                tf = dp*eta*xx-dq
!
! - XIF (T+DT)
!
!                CALL LCPRSC ( DFDS    , EPXI  , ZZ   )
                xx = zz*(1.d0-eta)*dp*3.d0/2.d0
                vtmp(1:ndt) = xx*epxi(1:ndt)
                xif(1:ndt) = vtmp(1:ndt)-dxxi(1:ndt)
!
            end if
        end if
!
! - RES (T+DT) = ( GF LF JF KF RF TF XIF (FF) )
!
!        PRINT *," TF ",TF
!        PRINT *," XIF ",XIF
!
        res(3*ndt+3) = tf
        res(3*ndt+4:3*ndt+4-1+ndt) = xif(1:ndt)
!
    end if
!
!
end subroutine
