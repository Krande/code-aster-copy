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
subroutine cvmini(typess, essai, mod, nmat, materf, &
                  timed, timef, yd, epsd, deps, &
                  dy)
    implicit none
!     VISCOCHABOCHE : CALCUL SOLUTION ESSAI ( SI IOPTIO = 0 )
!                       DY = ( DSIG DX1 DX2 DP DR 0  (DEPS3))
!          AVEC         Y =  ( SIG  X1  X2  P  R  Q  (EPS3) )
!     VISCOCHABOCHE : CALCUL SOLUTION ESSAI ( SI IOPTIO = 2 )
!                       DY = ( DSIG DX1 DX2 DP DR DQ DXXI (DEPS3))
!          AVEC         Y =  ( SIG  X1  X2  P  R  Q  XXI  (EPS3) )
!       ----------------------------------------------------------------
!       IN  ESSAI  :  VALEUR DE LA SOLUTION D ESSAI
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           TIMED  :  TEMPS A T
!           TIMEF  :  TEMPS A T+DT
!           YD     :  VARIABLES A T   =  ( SIG X1 X2 P R (EPS3) )
!                     OU (SI IOPTIO =2) ( SIG X1 X2 P R Q XXI (EPS3) )
!           EPSD   :  DEFORMATION A T
!       VAR DEPS   :  INCREMENT DE DEFORMATION
!           TYPESS :  TYPE DE SOLUTION D ESSAI
!                               0 = NUL(0)
!                               1 = ELASTIQUE
!                               2 = EXPLICITE (=-1 INITIALEMENT)
!                               3 = ESSAI
!       OUT DY     :  SOLUTION ESSAI  = ( DSIG DX1 DX2 DP DR (DEPS3))
!                            OU  ( SIG  X1  X2  P  R  Q  XXI  (EPS3))
!       ----------------------------------------------------------------
!
#include "asterfort/chbfs.h"
#include "asterfort/cvmcvx.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
    integer(kind=8) :: ndt, ndi, typess, nmat
    integer(kind=8) :: ioptio, idnr, nopt
!
!
    real(kind=8) :: yd(*), dy(*), essai
    real(kind=8) :: hook(6, 6), dfds(6), fkooh(6, 6)
    real(kind=8) :: epsd(6), deps(6)
    real(kind=8) :: sig(6), dsig(6)
    real(kind=8) :: x1(6), dx1(6), x2(6), dx2(6)
    real(kind=8) :: xxi(6), dxxi(6), p, dp
    real(kind=8) :: r, dr, q, dq
!
    real(kind=8) :: timed, timef, dt, seuil
    real(kind=8) :: materf(nmat, 2)
!
    real(kind=8) :: k0, ak, n, alp
    real(kind=8) :: b, mr, gr, mu, qm, q0
    real(kind=8) :: qr0, eta, ai
    real(kind=8) :: m1, d1, gx1, g10, c1, c1d
    real(kind=8) :: m2, d2, gx2, g20, c2, c2d
    real(kind=8) :: ccin, xx, yy, zz, nun, sgn
    real(kind=8) :: grq, qr
    real(kind=8) :: vtmp(6), vtmp1(6), epsp(6), epxi(6)
!
!
    character(len=8) :: mod
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/opti/ioptio, idnr
    common/coed/c1d, c2d
!       ----------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8) :: difc1, difc2
!-----------------------------------------------------------------------
    if (typess .eq. -1) typess = 2
!
!
    sig(1:ndt) = yd(1:ndt)
    x1(1:ndt) = yd(ndt+1:ndt+ndt)
    x2(1:ndt) = yd(2*ndt+1:2*ndt+ndt)
    p = yd(3*ndt+1)
    r = yd(3*ndt+2)
    q = yd(3*ndt+3)
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
    nopt = 0
    if (ioptio .eq. 2) nopt = idnr
!
    call lcopli('ISOTROPE', mod, materf(1, 1), hook)
!
! - SOLUTION INITIALE = NUL
!
    if (typess .eq. 0) then
        dy(1:4*ndt+4) = 0.d0
!
! - SOLUTION INITIALE = ELASTIQUE
!
    else if (typess .eq. 1) then
        dy(1:4*ndt+4) = 0.d0
        dsig(1:ndt) = matmul(hook(1:ndt, 1:ndt), deps(1:ndt))
        dy(1:ndt) = dsig(1:ndt)
!
! - SOLUTION INITIALE = EXPLICITE
!
    else if (typess .eq. 2) then
! -     SOLUTION D'ESSAI POUR  ( SIG  X1  X2  P  R  0  (EPS3))
!
        ccin = ai+(1.d0-ai)*exp(-b*p)
        dt = timef-timed
! - DP
        call cvmcvx(nmat, materf, sig, yd(ndt+1), seuil)
        if (seuil .le. 0.d0) then
            dp = 0.d0
            do i = 1, 6
                dfds(i) = 0.d0
            end do
        else
            call chbfs(sig, x1, x2, dfds)
            zz = seuil/(k0+ak*r)
!              IF ( ZZ .LT. 0.D0 ) ZZ = 0.D0
            dp = dt*(zz**n)*exp(alp*(zz**(n+1)))
        end if
!
! - DSIG
        vtmp(1:ndt) = (-dp)*dfds(1:ndt)
        vtmp(1:ndt) = deps(1:ndt)+vtmp(1:ndt)
        dsig(1:ndt) = matmul(hook(1:ndt, 1:ndt), vtmp(1:ndt))
!
! - DX1
        zz = dot_product(x1(1:ndt), dfds(1:ndt))
        yy = dot_product(x1(1:ndt), x1(1:ndt))
        zz = zz*(1.d0-d1)*g10*ccin*dp*2.d0/3.d0
        xx = c1*dp*2.d0/3.d0-zz
        if (yy .le. 0.d0) then
            yy = 0.d0
        else
            yy = gx1*dt*(sqrt(yy*3.d0/2.d0))**(m1-1.d0)+g10*ccin*d1*dp
        end if
        vtmp(1:ndt) = yy*x1(1:ndt)
        dx1(1:ndt) = xx*dfds(1:ndt)
        dx1(1:ndt) = dx1(1:ndt)-vtmp(1:ndt)
!
! - CAS ANISOTHERME
!
        if (c1 .ne. 0.d0) then
            difc1 = (c1-c1d)/c1
            vtmp(1:ndt) = difc1*x1(1:ndt)
            dx1(1:ndt) = dx1(1:ndt)+vtmp(1:ndt)
        end if
!
! - DX2
        zz = dot_product(x2(1:ndt), dfds(1:ndt))
        yy = dot_product(x2(1:ndt), x2(1:ndt))
        zz = zz*(1.d0-d2)*g20*ccin*dp*2.d0/3.d0
        xx = c2*dp*2.d0/3.d0-zz
        if (yy .le. 0.d0) then
            yy = 0.d0
        else
            yy = gx2*dt*(sqrt(yy*3.d0/2.d0))**(m2-1.d0)+g20*ccin*d2*dp
        end if
        vtmp(1:ndt) = yy*x2(1:ndt)
        dx2(1:ndt) = xx*dfds(1:ndt)
        dx2(1:ndt) = dx2(1:ndt)-vtmp(1:ndt)
!
! - CAS ANISOTHERME
!
        if (c2 .ne. 0.d0) then
            difc2 = (c2-c2d)/c2
            vtmp(1:ndt) = difc2*x2(1:ndt)
            dx2(1:ndt) = dx2(1:ndt)+vtmp(1:ndt)
        end if
!
! - DR
        grq = q0+(qm-q0)*(1.d0-exp(-2.d0*mu*q))
        qr = grq-qr0*(1.d0-((qm-grq)/qm)**2)
        sgn = 1.d0
        if ((qr-r) .lt. 0.d0) sgn = -1.d0
        dr = b*(grq-r)*dp+sgn*gr*dt*(abs(qr-r))**mr
!
! - DEPS(3)
        if (mod(1:6) .eq. 'C_PLAN') then
            nun = materf(2, 1)/(1.d0-materf(2, 1))
            deps(3) = nun*(dp*(dfds(1)+dfds(2))-deps(1)-deps(2))+dfds(3)*dp
        end if
!
!
! - DY
        dy(1:ndt) = dsig(1:ndt)
        dy(ndt+1:ndt+ndt) = dx1(1:ndt)
        dy(2*ndt+1:2*ndt+ndt) = dx2(1:ndt)
        dy(3*ndt+1) = dp
        dy(3*ndt+2) = dr
        dy(3*ndt+3) = 0.d0
        if (mod(1:6) .eq. 'C_PLAN') then
            dy(3*ndt+4+nopt) = deps(3)
            dy(3) = 0.d0
        end if
!
!
! -         SOLUTION D'ESSAI POUR ( SIG  X1  X2  P  R  Q  XXI  (EPS3))
        if (ioptio .eq. 2) then
            xxi(1:ndt) = yd(3*ndt+4:3*ndt+4-1+ndt)
!
! - EPSP
            call lcopil('ISOTROPE', mod, materf(1, 1), fkooh)
            vtmp(1:ndt) = sig(1:ndt)+dsig(1:ndt)
            vtmp1(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), vtmp(1:ndt))
            epsp(1:ndt) = epsd(1:ndt)+deps(1:ndt)
            epsp(1:ndt) = epsp(1:ndt)-vtmp1(1:ndt)
!
! N-ETOILE
            vtmp(1:ndt) = epsp(1:ndt)-xxi(1:ndt)
            xx = dot_product(vtmp(1:ndt), vtmp(1:ndt))
            xx = sqrt(xx*3.d0/2.d0)
!
! H(F)
            zz = 2.d0/3.d0*xx-q
!
            if (zz .lt. 0.d0) then
                dq = 0.d0
                dxxi(:) = 0.d0
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
                    dq = 0.d0
                    dxxi(:) = 0.d0
                else
!
! - DXXI
!                   CALL LCPRSC ( DFDS    , EPXI  , ZZ    )
                    xx = zz*(1.d0-eta)*dp*3.d0/2.d0
                    dxxi(1:ndt) = xx*epxi(1:ndt)
!
! - DQ
                    xx = dot_product(dfds(1:ndt), epxi(1:ndt))
                    dq = dp*eta*xx
                end if
            end if
!
! - DY
            dy(3*ndt+3) = dq
            dy(3*ndt+4:3*ndt+4-1+ndt) = dxxi(1:ndt)
        end if
!
! - SOLUTION INITIALE = VALEUR ESSAI POUR TOUTES LES COMPOSANTES
!
    else if (typess .eq. 3) then
        dy(1:4*ndt+4) = essai
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = essai
            dy(3) = 0.d0
        end if
!
    end if
!
!
end subroutine
