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

subroutine rsljpl(fami, kpg, ksp, loi, imat, &
                  nmat, mater, sig, vin, vind, &
                  deps, theta, dt, dsde)
    implicit none
!       ROUSSELIER :  MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT
!            COHERENT ELASTO-PLASTIQUE EN VITESSE A T+DT OU T
!       IN  FAMI   :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!           KPG,KSP:  NUMERO DU (SOUS)POINT DE GAUSS
!           IMAT   :  ADRESSE DU MATERIAU CODE
!           NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!           SIG    :  CONTRAINTES
!           VIN    :  VARIABLES INTERNES
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!       ----------------------------------------------------------------
#include "asterfort/lchydr.h"
#include "asterfort/lcnrte.h"
#include "asterfort/lcnrts.h"
#include "asterfort/lcprte.h"
#include "asterfort/lcsomh.h"
#include "asterfort/rsliso.h"
    integer(kind=8) :: kpg, ksp, imat, nmat
!
    real(kind=8) :: v1(6), v2(6), i2(6)
    real(kind=8) :: m1(6, 6), m2(6, 6), m3(6, 6), dsde(6, 6), i4(6, 6)
    real(kind=8) :: vin(*), vind(*)
    real(kind=8) :: deps(6), sig(6), rig(6), rigdv(6)
    real(kind=8) :: rigmo, rigeq, rigeq2, eps0, sig0, mexpo, acc
    real(kind=8) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, expo, expe
    real(kind=8) :: x1, x2, y1, y2, y3, y4, y5, a1, a2, a3, a4
    real(kind=8) :: a5, a6, d, s1, f, fd, f0, p, rho, unf
    real(kind=8) :: nu, e, deuxmu, mu, troimu, troisk, k
    real(kind=8) :: mater(nmat, 2), rp, drdp, dp, dpm
    real(kind=8) :: zero, un, deux, trois, ann, pd, ftot, fdtot
    real(kind=8) :: ndeps, nsig, theta
    real(kind=8) :: dpuiss, puiss, dt
!
    character(len=16) :: loi
    character(len=*) :: fami
!
    integer(kind=8) :: ndt, ndi
    common/tdim/ndt, ndi
!
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
!
    data i4/un, zero, zero, zero, zero, zero,&
     &                   zero, un, zero, zero, zero, zero,&
     &                   zero, zero, un, zero, zero, zero,&
     &                   zero, zero, zero, un, zero, zero,&
     &                   zero, zero, zero, zero, un, zero,&
     &                   zero, zero, zero, zero, zero, un/
    data i2/un, un, un, zero, zero, zero/
!       ----------------------------------------------------------------
!
! -- INITIALISATION----------------------------------------------
!
    p = vin(1)
    f = vin(2)
    pd = vind(1)
    fd = vind(2)
    e = mater(1, 1)
    nu = mater(2, 1)
    d = mater(1, 2)
    s1 = mater(2, 2)
    f0 = mater(3, 2)
    if (loi(1:10) .eq. 'ROUSS_VISC') then
        ann = 0.d0
        sig0 = mater(9, 2)
        eps0 = mater(10, 2)
        mexpo = mater(11, 2)
    else if (loi(1:10) .eq. 'ROUSS_PR') then
        ann = mater(8, 2)
        sig0 = 0.d0
        eps0 = 0.d0
        mexpo = 0.d0
    end if
    ftot = f+ann*p
    fdtot = fd+ann*pd
!
!
! -- CAS DU MATERIAU CASSE---------------------------------------
    if (fdtot .ge. mater(6, 2)) then
        ndeps = lcnrte(deps)
        nsig = lcnrts(sig)
        if ((ndeps*nsig) .eq. zero) then
            dsde(:, :) = zero
        else
            a1 = -deux*mater(7, 2)*e/(ndeps*nsig*trois)
            call lcprte(sig, deps, dsde)
            dsde(1:ndt, 1:ndt) = a1*dsde(1:ndt, 1:ndt)
        end if
!
! -- CAS DU MATERIAU NON CASSE-----------------------------------
    else
!
        call rsliso(fami, kpg, ksp, '-', imat, &
                    p, rp, drdp)
!
        unf = un-f
        rho = (unf-ann*p)/(un-f0)
        rig(1:ndt) = (un/rho)*sig(1:ndt)
        call lchydr(rig, rigmo)
        call lcsomh(rig, -rigmo, rigdv)
        rigeq = lcnrts(rigdv)
        rigeq2 = rigeq*rigeq
!
        deuxmu = e/(un+nu)
        mu = deuxmu/deux
        troimu = trois*mu
        troisk = e/(un-deux*nu)
        k = troisk/trois
!
! -- SI LA POROSITE EST ACCELERE----
        if (fdtot .ge. mater(4, 2)) then
            acc = mater(5, 2)
! -- SINON-----------------------
        else
            acc = un
        end if
! -- ----------------------------
!
        dp = p-vind(1)
        if (dp .gt. zero) then
            dpm = theta*(f-fd)/(trois*unf*acc)
            expo = d*exp(rigmo/s1)
            expe = expo*s1*ftot
!
            if (loi(1:10) .eq. 'ROUSS_VISC') then
                puiss = (dp/(dt*eps0))**(un/mexpo)
                dpuiss = ((dp/(dt*eps0))**(un/mexpo-un))/(mexpo*(dt*eps0))
                drdp = drdp+sig0*dpuiss/sqrt(un+puiss**2)/theta
            end if
!
            z1 = un+trois*dpm*acc
            z2 = troimu+drdp
            z3 = k*ftot*z1-s1*unf*acc
            z4 = dp*theta*drdp-rigeq
            z5 = troimu*theta*dp+rigeq
!
!
            x1 = expe*expo*z3+expo*z2*z3*theta*dp+z1*z2*s1-ann*z1*expo*s1*s1
            x2 = -(expe*expo*z3*dp+expo*z3*z4*theta*dp+z1*s1*z4)+ann*z1*expo*s1*s1*dp
!
            y1 = -troisk*expo*z1*ftot/(x1*rigeq)
            y2 = -troimu/(x1*z5*rigeq2)
            y3 = -troisk*expo*z1*ann*theta*dp/(x1*rigeq)
!
            a1 = troisk+y1*k*rigeq*(expe+dp*theta*z2)
            a2 = mu*(y1+y3)*s1
            a3 = deuxmu*(un+y2*x1*dp*theta*rigeq2)
            a4 = y2*troimu*x2
            a5 = y2*troisk*expe*z1*z5*rigeq
!
            z6 = expo
            z7 = z6*s1*ftot
            z8 = unf/(unf-ann*p)
            z9 = ann/(unf-ann*p)
            a6 = troimu*k*theta*dp-a2*rigeq*s1
            y4 = acc*z8/z1+z9*s1/(z7+z2*theta*dp)
            y5 = acc*a2*z8/z1-z9*a6/(rigeq*(z7+z2*theta*dp))
!
            m1(1:ndt, 1:ndt) = a3*i4(1:ndt, 1:ndt)
            v1(1:ndt) = ((a1-a3)/trois)*i2(1:ndt)
            v2(1:ndt) = a2*rigdv(1:ndt)
            v1(1:ndt) = v1(1:ndt)+v2(1:ndt)
            call lcprte(i2, v1, m2)
            v1(1:ndt) = a4*rigdv(1:ndt)
            v2(1:ndt) = (a5/trois)*i2(1:ndt)
            v1(1:ndt) = v1(1:ndt)+v2(1:ndt)
            call lcprte(rigdv, v1, m3)
            dsde(1:ndt, 1:ndt) = m1(1:ndt, 1:ndt)+m2(1:ndt, 1:ndt)+m3(1:ndt, 1:ndt)
!
! A CE STADE DSDE EST LE TENSEUR TANGENT COHERENT
! ENTRE D(SIG/RHO) ET DEPS
!
            v1(1:ndt) = (a1-troisk)*i2(1:ndt)
            v2(1:ndt) = (trois*y5/y4)*rigdv(1:ndt)
            v1(1:ndt) = v1(1:ndt)+v2(1:ndt)
            call lcprte(rig, v1, m1)
            m1(1:ndt, 1:ndt) = (y4/troisk)*m1(1:ndt, 1:ndt)
            dsde(1:ndt, 1:ndt) = dsde(1:ndt, 1:ndt)+m1(1:ndt, 1:ndt)
            dsde(1:ndt, 1:ndt) = rho*dsde(1:ndt, 1:ndt)
! -- CAS DP=0
        else
            call lcprte(i2, i2, m1)
            m1(1:ndt, 1:ndt) = (un/trois)*m1(1:ndt, 1:ndt)
            m2(1:ndt, 1:ndt) = (rho*troisk)*m1(1:ndt, 1:ndt)
            m1(1:ndt, 1:ndt) = (-un)*m1(1:ndt, 1:ndt)
            m1(1:ndt, 1:ndt) = i4(1:ndt, 1:ndt)+m1(1:ndt, 1:ndt)
            m1(1:ndt, 1:ndt) = (rho*deuxmu)*m1(1:ndt, 1:ndt)
            dsde(1:ndt, 1:ndt) = m1(1:ndt, 1:ndt)+m2(1:ndt, 1:ndt)
        end if
    end if
! ------------------------------------------------------------------
end subroutine
