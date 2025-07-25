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

subroutine cjstde(mod, mater, nvi, eps, sig, &
                  vin, dsde)
    implicit none
!     CALCUL DE LA MATRICE TANGENTE DU PROBLEME CONTINU DE LA LOI CJS
!     POUR LE MECANISME PLASTIQUE DEVIATOIRE
!     IN   MOD     :  MODELISATION
!          MATER   :  COEFFICIENTS MATERIAU
!          NVI     :  NB DE VARIABLES INTERNES
!          EPS     :  DEFORMATIONS
!          SIG     :  CONTRAINTES
!          VIN     :  VARIABLES INTERNES
!     OUT  DSDESY  :  MATRICE TANGENTE
!       ----------------------------------------------------------------
! ======================================================================
#include "asterfort/calcq.h"
#include "asterfort/cjsqij.h"
#include "asterfort/cos3t.h"
#include "asterfort/hlode.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lcdevi.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nvi, i, j, codret
    real(kind=8) :: mater(14, 2), vin(*), sig(6), dsde(6, 6)
    real(kind=8) :: hook(6, 6), i1, eps(6), epsv, e, nu, al, la, mu
    real(kind=8) :: r, gr, x(6), gx(6), xii, dfdds(6), gd(6), kron(6)
    real(kind=8) :: s(6), sii, hts, siic, norm(6), epssig, siirel, pref
    real(kind=8) :: n, rm, rc, a, b, c, pco, pa, qinit, hdev, qgx
    real(kind=8) :: beta, betapr, gamma, mucjs, dfhgd, dfh(6), hgd(6)
    real(kind=8) :: q(6), qii, qq(6), qqii
    real(kind=8) :: pc, cos3tq, htq, cosa, cosdif, rr, phi, phio
    real(kind=8) :: cos3ts, tangs, tangq, tetas, tetaq
    real(kind=8) :: truc, signe, coef3, coef4
    real(kind=8) :: mun5, zero, d12, un, deux, trois
    character(len=8) :: mod
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
! --- DEFINITION DE PARAMETRES -----------------------------------------
! ======================================================================
    parameter(epssig=1.d-8)
! ======================================================================
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data mun5/-1.5d0/
    data zero/0.d0/
    data d12/.5d0/
    data un/1.d0/
    data deux/2.d0/
    data trois/3.d0/
! ======================================================================
    call jemarq()
! ======================================================================
! --- RECUPERATION DES GRANDEURS UTILES : I1, VARIABLES INTERNES R ET X,
! --- PARAMETRES CJS, RIGIDITE HOOK ------------------------------------
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
    qinit = mater(13, 2)
! ======================================================================
! --- CALCUL DE LA TRACE DE SIG ----------------------------------------
! ======================================================================
    i1 = trace(ndi, sig)
    if ((i1+qinit) .eq. 0.0d0) then
        i1 = -qinit+1.d-12*pa
        pref = abs(pa)
    else
        pref = abs(i1+qinit)
    end if
! ======================================================================
! --- RECUPERATION DES VARIABLES INTERNES ------------------------------
! ======================================================================
    r = vin(2)
    do i = 1, ndt
        x(i) = vin(i+2)
    end do
! ======================================================================
! --- OPERATEUR DE RIGIDITE CALCULE A ITERARTION -----------------------
! ======================================================================
    e = mater(1, 1)*((i1+qinit)/trois/pa)**n
    nu = mater(2, 1)
    al = e*(un-nu)/(un+nu)/(un-deux*nu)
    la = nu*e/(un+nu)/(un-deux*nu)
    mu = e*d12/(un+nu)
    hook(:, :) = zero
! ======================================================================
! --- 3D/DP/AX ---------------------------------------------------------
! ======================================================================
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
        ! ======================================================================
        ! --- CP/1D ------------------------------------------------------------
        ! ======================================================================
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'ALGORITH2_15')
        ! ======================================================================
    end if
! ======================================================================
! --- LOIS D'ECOUISSAGE DES VARIABLES INTERNES R ET X ------------------
! ======================================================================
! --- ECROUISSAGE ISOTROPE ---------------------------------------------
! ======================================================================
    gr = -a*(un-r/rm)**deux*(i1+qinit)*((i1+qinit)/trois/pa)**mun5
! ======================================================================
! --- ECROUISSAGE CINEMATIQUE ------------------------------------------
! ======================================================================
    call lcdevi(sig, s)
    sii = norm2(s(1:ndt))
    siirel = sii/pref
    cos3ts = cos3t(s, pref, epssig)
    hts = hlode(gamma, cos3ts)
!
    call cjsqij(s, i1, x, q)
    qii = norm2(q(1:ndt))
    cos3tq = cos3t(q, pref, epssig)
    htq = hlode(gamma, cos3tq)
! ======================================================================
! --- CALCUL DE Q ------------------------------------------------------
! ======================================================================
    call calcq(q, gamma, pref, epssig, qq, &
               codret)
! ======================================================================
! --- CALCUL DE PC (CONTRAINTE MOYENNE CRITIQUE) -----------------------
! ======================================================================
    qqii = norm2(qq(1:ndt))
    xii = norm2(x(1:ndt))
    epsv = trace(ndi, eps)
    pc = pco*exp(-c*epsv)
! ======================================================================
! --- CALCUL DE LA FONCTION PHI ----------------------------------------
! ======================================================================
    if (xii .le. epssig) then
        phi = un
    else if (siirel .le. epssig) then
        cosa = un
        cosdif = un
        rr = rc+mucjs*max(zero, log(trois*pc/(i1+qinit)))
        phio = cosa/(rr-hts/htq*rm*cosdif)
        phi = phio*hts*qqii
    else
        cosa = (qii*qii-sii*sii-i1*xii*i1*xii)/(deux*sii*i1*xii)
        tangs = sqrt(un-cos3ts*cos3ts)/cos3ts
        tangq = sqrt(un-cos3tq*cos3tq)/cos3tq
        tetas = atan2(tangs, 1.d0)/trois
        tetaq = atan2(tangq, 1.d0)/trois
        cosdif = cos(tetas-tetaq)
        rr = rc+mucjs*max(zero, log(trois*pc/(i1+qinit)))
        phio = cosa/(rr-hts/htq*rm*cosdif)
        phi = phio*hts*qqii
    end if
! ======================================================================
! --- CALCUL DE GX -----------------------------------------------------
! ======================================================================
    do i = 1, ndt
        gx(i) = (i1+qinit)/b*(qq(i)+phi*x(i))*((i1+qinit)/trois/pa)**mun5
    end do
! ======================================================================
! --- LOI D'ECOULEMENT DU MECANISME DEVIATOIRE -------------------------
! ======================================================================
    truc = dot_product(qq(1:ndt), x(1:ndt))-r
    do i = 1, ndi
        dfdds(i) = qq(i)-truc
    end do
    do i = ndi+1, ndt
        dfdds(i) = qq(i)
    end do
    siic = -rc*(i1+qinit)/hts
    signe = vin(nvi-1)
    betapr = beta*(sii/siic-un)*signe
    coef4 = un/sqrt(betapr*betapr+trois)
    coef3 = coef4*betapr/sii
    do i = 1, ndt
        norm(i) = coef3*s(i)+coef4*kron(i)
    end do
    truc = dot_product(dfdds(1:ndt), norm(1:ndt))
    do i = 1, ndt
        gd(i) = dfdds(i)-truc*norm(i)
    end do
! ======================================================================
! --- MODULE PLASTIQUE DEVIATOIRE : HDEV -------------------------------
! ======================================================================
    qgx = zero
    do i = 1, ndt
        qgx = qgx+qq(i)*gx(i)
    end do
    hdev = (i1+qinit)*(-gr+qgx)
! ======================================================================
! --- CALCUL DU TERME DFDDS.HOOK.GD : DFHGD ----------------------------
! ======================================================================
    dfhgd = zero
    do i = 1, ndt
        do j = 1, ndt
            dfhgd = dfhgd+dfdds(i)*hook(i, j)*gd(j)
        end do
    end do

! ======================================================================
! --- CALCUL DU TERME HOOK.GD : HGD ------------------------------------
! ======================================================================
    hgd(:) = zero
    do i = 1, ndt
        do j = 1, ndt
            hgd(i) = hgd(i)+hook(i, j)*gd(j)
        end do

    end do
! ======================================================================
! --- CALCUL DU TERME DFDDS.HOOK : DFH ---------------------------------
! ======================================================================
    dfh(:) = zero
    do i = 1, ndt
        do j = 1, ndt
            dfh(i) = dfh(i)+dfdds(j)*hook(j, i)
        end do

    end do
! ======================================================================
! --- CALCUL DE DSDE(I,J,K,L) =
!
!                      HOOK(I,J,M,N) GD(M,N) DFDDS(P,Q) HOOK(P,Q,K,L)
!    HOOK(I,J,K,L) -  -----------------------------------------------
!                         HDEV  +  DFDDS(R,S) HOOK(R,S,T,U) GD(T,U)
!
! --- C'EST A DIRE :  DSDE = HOOK - HGD.DFH / ( HDEV + DFHGD ) ---------
! ======================================================================
    dsde(:, :) = 0.d0
    do i = 1, ndt
        do j = 1, ndt
            dsde(i, j) = hook(i, j)-hgd(i)*dfh(j)/(hdev+dfhgd)
        end do
    end do
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
