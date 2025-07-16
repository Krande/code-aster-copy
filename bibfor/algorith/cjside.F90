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
subroutine cjside(mod, mater, epsd, deps, yd, &
                  gd, dy)
    implicit none
!     LOI CJS :  MECANISME DEVIATOIRE
!     CALCUL DE LA SOLUTION D ESSAI EXPLICITE DY = (DSIG, DVIN, DLAMB )
!                           AVEC               Y = ( SIG,  VIN,  LAMB )
!     PARTIR DE LA PREDICTION ELASTIQUE
!     ------------------------------------------------------------------
!     IN   MOD      :  MODELISATION
!          MATER    :  COEFFICIENTS MATERIAU A T+DT
!          EPSD     :  DEFORMATION A T+DT
!          DEPS     :  INCREMENT DE DEFORMATION
!          YD       :  VARIABLES A T = (SIGD, VIND, LAMB)
!     VAR  GD       :  TENSEUR DE LA LOI D ECOULEMENT PLASTIQUE DEV.
!     OUT  DY       :  SOLUTION D ESSAI (DSIG, DVIN, DLAMB)
!     ------------------------------------------------------------------
!     Y CONTIENT LES CONTRAINTES : SIG
!                LES VARIABLES INTERNES :  R, X
!                LE MULTIPLICATEUR PLASTIQUE :  LAMBD
! ======================================================================
#include "asterfort/calcq.h"
#include "asterfort/cjsqco.h"
#include "asterfort/cjst.h"
#include "asterfort/cos3t.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, i, j, codret
    real(kind=8) :: epsd(6), deps(6), depse(6), trdeps, kron(6)
    real(kind=8) :: dsig(6), sigd(6), sige(6), yd(*), dy(*), dqe(6)
    real(kind=8) :: mater(14, 2), n, rm, rc, a, b, c, pco, pc, pa
    real(kind=8) :: ke, koe, beta, betapr, gamma, mucjs, dq(6)
    real(kind=8) :: hooknl(6, 6), epssig, pref, qinit
    real(kind=8) :: e, nu, al, la, mu, i1d, i1e
    real(kind=8) :: truc, signe, dlambd, rcos3t
    real(kind=8) :: s(6), sii, siic, hts, dets, cos3ts, siirel
    real(kind=8) :: se(6), siie, htse, detse, co3tse, siiere
    real(kind=8) :: q(6), qii, htq, detq, cos3tq, qiirel
    real(kind=8) :: qe(6), qiie, htqe, detqe, co3tqe, qiiere
    real(kind=8) :: tangs, tangq, tetas, tetaq
    real(kind=8) :: qq(6), qqii, norm(6), dfdds(6), gd(6), trgd
    real(kind=8) :: coef1, coef3, coef4, coef5, coef6, prod1
    real(kind=8) :: rd, xd(6), gr, gx(6), xii, dr, dx(6)
    real(kind=8) :: epsv, phi, phio, rr, cosa, cosdif
    real(kind=8) :: fd, dfddl, drdl, di1dl, dqiidl, dhdl
    real(kind=8) :: dqdl(6), mun5, zero, un, d12, deux, trois, cinq
    character(len=8) :: mod
! ======================================================================
    parameter(mun5=-1.5d0)
    parameter(d12=.5d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    parameter(cinq=5.d0)
    parameter(epssig=1.d-8)
! ======================================================================
    common/tdim/ndt, ndi
    data kron/un, un, un, zero, zero, zero/
! ======================================================================
    call jemarq()
! ======================================================================
! --- PROPRIETES CJS MATERIAU ------------------------------------------
! ======================================================================
    beta = mater(1, 2)
    rm = mater(2, 2)
    n = mater(3, 2)
    rc = mater(5, 2)
    a = mater(6, 2)
    b = mater(7, 2)
    c = mater(8, 2)
    gamma = mater(9, 2)
    mucjs = mater(10, 2)
    pco = mater(11, 2)
    pa = mater(12, 2)
    koe = mater(1, 1)/trois/(un-deux*mater(2, 1))
    qinit = mater(13, 2)
! ======================================================================
! --- PREMIER INVARIANT ET AUTRES GRANDEURS UTILES ---------------------
! ======================================================================
    i1d = trace(ndi, yd)
    if ((i1d+qinit) .eq. 0.0d0) then
        i1d = -qinit+1.d-12*pa
        pref = abs(pa)
    else
        pref = abs(i1d+qinit)
    end if
!
    rd = yd(ndt+1)
    do i = 1, ndt
        xd(i) = yd(ndt+1+i)
    end do
    ke = koe*((i1d+qinit)/trois/pa)**n
! ======================================================================
! --- OPERATEUR DE RIGIDITE NON LINEAIRE -------------------------------
! ======================================================================
! --- OPERATEUR LINEAIRE NON LINEAIRE ----------------------------------
! ======================================================================
    hooknl(:, :) = zero
    e = mater(1, 1)*((i1d+qinit)/trois/pa)**n
    nu = mater(2, 1)
    al = e*(un-nu)/(un+nu)/(un-deux*nu)
    la = nu*e/(un+nu)/(un-deux*nu)
    mu = e*d12/(un+nu)
! ======================================================================
! --- 3D/DP/AX ---------------------------------------------------------
! ======================================================================
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        do i = 1, ndi
            do j = 1, ndi
                if (i .eq. j) hooknl(i, j) = al
                if (i .ne. j) hooknl(i, j) = la
            end do
        end do
        do i = ndi+1, ndt
            do j = ndi+1, ndt
                if (i .eq. j) hooknl(i, j) = deux*mu
            end do
        end do
! ======================================================================
! --- CP/1D ------------------------------------------------------------
! ======================================================================
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'ALGORITH2_15')
    end if
! ======================================================================
! --- LOIS D ECROUISSAGE : GR ET GX ------------------------------------
! ======================================================================
! --- ECROUISSAGE ISOTROPE ---------------------------------------------
! ======================================================================
    coef1 = ((i1d+qinit)/trois/pa)**mun5
    gr = -a*(un-rd/rm)**deux*(i1d+qinit)*coef1
! ======================================================================
! --- ECROUISSAGE CINEMATIQUE ------------------------------------------
! ======================================================================
    do i = 1, ndt
        sigd(i) = yd(i)
    end do
! ======================================================================
! --- CALCUL DE S, SII, COS3TS, .... -----------------------------------
! ======================================================================
    call cjsqco(gamma, sigd, xd, pref, epssig, &
                i1d, s, sii, siirel, cos3ts, &
                hts, dets, q, qii, qiirel, &
                cos3tq, htq, detq)
! ======================================================================
! --- ON CALCULE DE TOUTES FACONS UNE PREDICTION ELASTIQUE -------------
! ======================================================================
    dsig(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), deps(1:ndt))
    sige(1:ndt) = sigd(1:ndt)+dsig(1:ndt)
    i1e = trace(ndi, sige)
    if ((i1e+qinit) .eq. zero) then
        i1e = -qinit+1.d-12*pa
    end if
!
    call cjsqco(gamma, sige, xd, pref, epssig, &
                i1e, se, siie, siiere, co3tse, &
                htse, detse, qe, qiie, qiiere, &
                co3tqe, htqe, detqe)
!
! ======================================================================
! --- SI QII EST QUASI-NULL, IL N'Y A PAS DE DEVIATEUR. ----------------
! --- DONC LE TENSEUR QQ(SIGD) N'EXISTE PAS. ON PRENDRA ALORS ----------
! --- A LA PLACE QQ(SIG_PREDICTION ELAS) -------------------------------
! ======================================================================
    if (qiirel .le. epssig) then
        call calcq(qe, gamma, pref, epssig, qq, &
                   codret)
    else
        call calcq(q, gamma, pref, epssig, qq, &
                   codret)
    end if
    qqii = norm2(qq(1:ndt))
!
    xii = norm2(xd(1:ndt))
!
    epsv = zero
    do i = 1, ndi
        epsv = epsv+epsd(i)+deps(i)
    end do
!
    pc = pco*exp(-c*epsv)
!
    if (xii .le. epssig) then
        phi = un
    else if (siirel .le. epssig) then
        cosa = un
        cosdif = un
        rr = rc+mucjs*max(zero, log(trois*pc/(i1d+qinit)))
        phio = cosa/(rr-hts/htq*rm*cosdif)
        phi = phio*hts*qqii
    else
        cosa = (qii*qii-sii*sii-i1d*i1d*xii*xii)/(deux*sii*i1d*xii)
!
        tangs = sqrt(un-cos3ts*cos3ts)/cos3ts
        tangq = sqrt(un-cos3tq*cos3tq)/cos3tq
        tetas = atan2(tangs, 1.d0)/trois
        tetaq = atan2(tangq, 1.d0)/trois
        cosdif = cos(tetas-tetaq)
!
        rr = rc+mucjs*max(zero, log(trois*pc/(i1d+qinit)))
        phio = cosa/(rr-hts/htq*rm*cosdif)
        phi = phio*hts*qqii
    end if
!
    do i = 1, ndt
        gx(i) = (i1d+qinit)/b*(qq(i)+phi*xd(i))*coef1
    end do
! ======================================================================
! --- LOI D ECOULEMENT : GD --------------------------------------------
! ======================================================================
    truc = dot_product(qq(1:ndt), xd(1:ndt))-rd
!
    do i = 1, ndt
        dfdds(i) = qq(i)-truc*kron(i)
    end do
! ======================================================================
! --- HYPOTHESE : SIGNE(S,DEPS) = SIGNE(S,DEPSDP) ----------------------
! ======================================================================
    truc = dot_product(s(1:ndt), deps(1:ndt))
    if (truc .ge. zero) then
        signe = un
    else
        signe = -un
    end if
!
    siic = -rc*(i1d+qinit)/hts
    betapr = beta*(sii/siic-un)*signe
!
    if (siirel .gt. epssig) then
        coef3 = betapr/sii
        coef4 = un/sqrt(betapr*betapr+trois)
!
        do i = 1, ndt
            norm(i) = coef4*(coef3*s(i)+kron(i))
        end do
    else
        coef3 = betapr/siie
        coef4 = un/sqrt(betapr*betapr+trois)
        do i = 1, ndt
            norm(i) = coef4*(coef3*se(i)+kron(i))
        end do
    end if
!
    prod1 = dot_product(dfdds(1:ndt), norm(1:ndt))
    do i = 1, ndt
        gd(i) = dfdds(i)-prod1*norm(i)
    end do
! ======================================================================
! --- CALCUL DE DLAMBD -------------------------------------------------
! ======================================================================
    fd = qiie*htqe+rd*(i1e+qinit)
    trgd = trace(ndi, gd)
    trdeps = trace(ndi, deps)
    di1dl = -trois*ke*trgd
    drdl = gr
    dqdl(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), gd(1:ndt))
    do i = 1, ndt
        dqdl(i) = -dqdl(i)+ke*trgd*(kron(i)+trois*xd(i))-gx(i)*(i1d+trois*ke*trdeps&
                  & )
    end do
! ======================================================================
! --- POUR LE CALCUL DE DQIIDL, DHDL: ----------------------------------
! --- SI QII EST QUASI-NULL, IL N'Y A PAS DE DEVIATEUR. ----------------
! ======================================================================
    if (qiirel .le. epssig) then
        dqiidl = zero
        do i = 1, ndt
            dqiidl = dqiidl+qe(i)/qiie*dqdl(i)
        end do
!
        rcos3t = cos3t(qe, pref, epssig)
        call cjst(qe, dqe)
        coef5 = sqrt(trois/deux)*gamma/htqe**cinq/qiie**trois
        coef6 = -gamma*rcos3t/(deux*htqe**cinq*qiie**deux)
!
        dhdl = zero
        do i = 1, ndt
            dhdl = dhdl+(coef5*dqe(i)+coef6*qe(i))*dqdl(i)
        end do
        dfddl = htqe*dqiidl+qii*dhdl+rd*di1dl+(i1d+qinit)*drdl
    else
        dqiidl = zero
        do i = 1, ndt
            dqiidl = dqiidl+q(i)/qii*dqdl(i)
        end do
!
        rcos3t = cos3t(q, pref, epssig)
        call cjst(q, dq)
        coef5 = sqrt(trois/deux)*gamma/htq**cinq/qii**trois
        coef6 = -gamma*rcos3t/(deux*htq**cinq*qii**deux)
!
        dhdl = zero
        do i = 1, ndt
            dhdl = dhdl+(coef5*dq(i)+coef6*q(i))*dqdl(i)
        end do
        dfddl = htq*dqiidl+qii*dhdl+rd*di1dl+(i1d+qinit)*drdl
    end if
!
    dlambd = -fd/dfddl
! ======================================================================
! --- CALCUL DES INCREMENTS DE DEFORMATIONS ELASTIQUE ------------------
! ======================================================================
    do i = 1, ndt
        depse(i) = deps(i)-dlambd*gd(i)
    end do
! ======================================================================
! --- CALCUL INCREMENT DE CONTRAINTES  DSIG = HOOKNL.DEPSE -------------
! ======================================================================
    dsig(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), depse(1:ndt))
! ======================================================================
! --- CALCUL INCREMENT DE LA VARIABLE INTERNE R ------------------------
! ======================================================================
    dr = dlambd*drdl
! ======================================================================
! --- CALCUL INCREMENT DE LA VARIABLE INTERNE X ------------------------
! ======================================================================
    do i = 1, ndt
        dx(i) = dlambd*gx(i)
    end do
! ======================================================================
! --- SOLUTION D ESSAI -------------------------------------------------
! ======================================================================
    do i = 1, ndt
        dy(i) = dsig(i)
    end do
    dy(ndt+1) = dr
    do i = 1, ndt
        dy(ndt+1+i) = dx(i)
    end do
    dy(2*ndt+2) = dlambd
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
