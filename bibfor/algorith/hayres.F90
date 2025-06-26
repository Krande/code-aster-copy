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
subroutine hayres(mod, nmat, materd, materf, timed, &
                  timef, yd, yf, deps, dy, &
                  res, crit, iret)
    implicit none
!     HAYHURST :
!            CALCUL DES TERMES DU SYSTEME NL A RESOUDRE = RES(DY)
!                   DY  = ( DEPSILON_EL DP DH1 DH2 DPHI DD )
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
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/fgequi.h"
#include "asterfort/lcopli.h"
#include "blas/dscal.h"
    character(len=8) :: mod
    integer(kind=8) :: iret, itens, ndi, nmat, ndt, ndim
    real(kind=8) :: hookf(6, 6), res(10), crit(*), theta, alphad
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), timed, timef, deps(6), dt
    real(kind=8) :: dtot
    real(kind=8) :: yd(*), yf(*), dy(*), sigf(6), gh, dmg, dmgmi
    real(kind=8) :: depsp(6), devcum, decrou(2), ddmg, epsef(6), depsel(6), m13
    real(kind=8) :: ze, td, sinn, grj0, gh1, gh2, equi(17), rmin, sequid
    real(kind=8) :: eps0, pk, ph1, ph2, delta1, delta2, h1st, h2st, pkc, sig0
    real(kind=8) :: biga
    real(kind=8) :: trsig, grj2v, grj1, epsi, terme1, shmax, sequi, dddmg, dh1
    real(kind=8) :: dh2, dp
    blas_int :: b_incx, b_n
!     ----------------------------------------------------------------
    parameter(ze=0.0d0)
    parameter(td=1.5d0)
!
    common/tdim/ndt, ndi
!-----------------------------------------------------------------------
    theta = crit(4)
    gh1 = yd(8)
    gh2 = yd(9)
    dp = dy(7)
    dh1 = dy(8)
    dh2 = dy(9)
    dddmg = dy(10)
    do itens = 1, ndt
        epsef(itens) = yd(itens)+theta*dy(itens)
    end do
    dt = timef-timed
    if (ndt .eq. 6) then
        ndim = 3
    else if (ndt .eq. 4) then
        ndim = 2
        sigf(5) = 0.d0
        sigf(6) = 0.d0
        depsp(5) = 0.d0
        depsp(6) = 0.d0
    end if
    iret = 0
    rmin = r8miem()
    shmax = 50.d0
    eps0 = materf(1, 2)
    pk = materf(2, 2)
    ph1 = materf(3, 2)
    ph2 = materf(4, 2)
    delta1 = materf(5, 2)
    delta2 = materf(6, 2)
    h1st = materf(7, 2)
    h2st = materf(8, 2)
    biga = materf(9, 2)
    sig0 = materf(10, 2)
    pkc = materf(11, 2)
    alphad = materf(12, 2)
    sequid = materf(13, 2)
    epsi = r8prem()*pk
    gh = gh1+theta*dh1+gh2+theta*dh2
    m13 = -1.d0/3.d0
    dmgmi = 1.d0-(1.d0+pkc*timef)**m13
    dmg = yd(10)+theta*dddmg
!
!----------------------------------------------------------------
    dtot = (1.d0-dmg)
    call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
    sigf(1:ndt) = matmul(hookf(1:ndt, 1:ndt), epsef(1:ndt))
    sigf(1:ndt) = dtot*sigf(1:ndt)
!
!------------CALCUL DES INVARIANTS DE CONTRAINTE  -------
!     attention FGEQUI ne prend pas en compte les SQRT(2)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, 1.d0/sqrt(2.d0), sigf(4), b_incx)
    call fgequi(sigf, 'SIGM_DIR', ndim, equi)
!     on retablit le tenseur
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, sqrt(2.d0), sigf(4), b_incx)
    trsig = equi(16)
    grj0 = max(equi(3), equi(4))
    grj0 = max(grj0, equi(5))
    grj1 = trsig
    grj2v = equi(1)
    if (sequid .eq. 0.d0) then
        sequi = grj0
    else
        sequi = grj1
    end if
!------------ CALCUL DU TENSEUR DEPSPATORIQUE DES CONTRAINTES ---
    do itens = 1, ndt
        if (itens .le. 3) sigf(itens) = sigf(itens)-grj1/3.d0
    end do
!
!----- EQUATION DONNANT LA DERIVEE DE LA DEF VISCO PLAST
!----- CUMULEE
!
    terme1 = (grj2v*(1-gh))/(pk*(1-dmgmi)*dtot)
    if (grj2v .le. epsi) then
        devcum = ze
    else if (abs(terme1) .lt. shmax) then
        devcum = eps0*(sinh(terme1))*dt
    else
        iret = 1
        goto 999
    end if
!
!----- EQUATION DONNANT LA DERIVEE DE GH
!
    if (grj2v .le. epsi) then
!       DIVISION PAR ZERO EVITEE
        decrou(1) = ze
        decrou(2) = ze
    else
        if (gh1 .le. (h1st-rmin)) then
            decrou(1) = (ph1/grj2v)*(h1st-(delta1*(gh1+theta*dh1)))*dp
            decrou(2) = (ph2/grj2v)*(h2st-(delta2*(gh2+theta*dh2)))*dp
        else
            iret = 1
            goto 999
        end if
    end if
!
!----- EQUATION DONNANT LA DERIVEE DE L ENDOMMAGEMENT
!
    if (sequi .ge. ze) then
        sinn = alphad*sequi+((1.d0-alphad)*grj2v)
    else
        sinn = (1.d0-alphad)*grj2v
    end if
    if ((sinn/sig0) .lt. shmax) then
        ddmg = biga*sinh(sinn/sig0)*dt
    else
        iret = 1
        goto 999
    end if
!
!------ EQUATION DONNANT LA DERIVEE DE LA DEF VISCO PLAST
!
    if (grj2v .le. epsi) then
        do itens = 1, ndt
            depsp(itens) = ze
        end do
    else
        do itens = 1, ndt
            depsp(itens) = td*dp*sigf(itens)/grj2v
        end do
    end if
    depsel(1:ndt) = deps(1:ndt)-depsp(1:ndt)
    res(1:ndt) = depsel(1:ndt)-dy(1:ndt)
    res(ndt+1) = devcum-dp
    res(ndt+2) = decrou(1)-dh1
    res(ndt+3) = decrou(2)-dh2
    res(ndt+4) = ddmg-dddmg
999 continue
end subroutine
