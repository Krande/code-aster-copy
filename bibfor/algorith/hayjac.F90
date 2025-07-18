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
subroutine hayjac(mod, nmat, coefel, coeft, timed, &
                  timef, yf, deps, nr, nvi, &
                  vind, vinf, yd, dy, crit, &
                  drdy, iret)
    implicit none
!     --------------------------------------------------------------
!     CALCUL DU JACOBIEN DE HAYHURST = DRDY(DY)
!     IN  MOD    :  TYPE DE MODELISATION
!         NMAT   :  DIMENSION MATER
!         MATERF :  COEFFICIENTS MATERIAU A T+DT
!         YF     :  VARIABLES A T + DT =  ( SIGF DLAMBDA XI_P XI_VP)
!         DEPS   :  INCREMENT DE DEFORMATION
!         TIMED  :  INSTANT  T
!         TIMEF  :  INSTANT  T+DT
!         NR     :  DIMENSION DECLAREE DRDY
!         NVI    :  NOMBRE DE VARIABLES INTERNES
!         VIND   :  VARIABLE INTERNES A T
!         VINF   :  VARIABLE INTERNES A T+DT
!         YD     :  VARIABLES A T  = ( SIGD  0 XI_P XI_VP) A T
!         DY     :  SOLUTION = ( DSIG  DLAMBDA  DXI_P DXI_VP )
!     OUT DRDY   :  JACOBIEN DU SYSTEME NON LINEAIRE
!         IRET   :  CODE RETOUR
!     --------------------------------------------------------------
#include "asterc/r8prem.h"
#include "asterfort/fgequi.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcprte.h"
#include "asterfort/r8inir.h"
#include "blas/dscal.h"
    integer(kind=8) :: nr, nmat, iret, nvi, i, j
    real(kind=8) :: deps(6), drdy(nr, nr), yf(nr), dy(nr), yd(nr)
    real(kind=8) :: coefel(nmat), coeft(nmat), crit(*)
    real(kind=8) :: timed, timef, vind(nvi), vinf(nvi)
    character(len=8) :: mod
!
    integer(kind=8) :: ndt, ndi, itens, ndim
    real(kind=8) :: dev(6), n(6)
    real(kind=8) :: kron(6), eps0, pk, ph1, ph2, delta1, delta2, h1st, h2st
    real(kind=8) :: biga, sig0, pkc, alphad, sequid, dt, theta, epsi, coef
    real(kind=8) :: h1, h2, h, d, dp, cosh1, cosh2, sinh1, sinh2, dmgmi, grj0
    real(kind=8) :: seq, seq0, sequi, shmax, sinn, terme1, trsig, trsig0, ddmg
    real(kind=8) :: sigf(6), equi(16), dm1, epsef(6), young, poiss, deuxmu, dmg
    real(kind=8) :: id(6, 6), nxn(6, 6), dfedee(6, 6), dh1
    real(kind=8) :: dh2
    real(kind=8) :: un, zero, d23, d13, dseqde(6), hookf(6, 6), troisk, gh1, gh2
    blas_int :: b_incx, b_n
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(d23=2.d0/3.d0)
    parameter(d13=-1.d0/3.d0)
!
!     --------------------------------------------------------------
    common/tdim/ndt, ndi
!     --------------------------------------------------------------
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data id/d23, d13, d13, zero, zero, zero,&
     &                  d13, d23, d13, zero, zero, zero,&
     &                  d13, d13, d23, zero, zero, zero,&
     &                  zero, zero, zero, un, zero, zero,&
     &                  zero, zero, zero, zero, un, zero,&
     &                  zero, zero, zero, zero, zero, un/
!
! --    COEFFICIENTS MATERIAU INELASTIQUES
!
    eps0 = coeft(1)
    pk = coeft(2)
    ph1 = coeft(3)
    ph2 = coeft(4)
    delta1 = coeft(5)
    delta2 = coeft(6)
    h1st = coeft(7)
    h2st = coeft(8)
    biga = coeft(9)
    sig0 = coeft(10)
    pkc = coeft(11)
    alphad = coeft(12)
    sequid = coeft(13)
    young = coefel(1)
    poiss = coefel(2)
    deuxmu = young/(1.d0+poiss)
    troisk = young/(1.d0-2.d0*poiss)
    dt = timef-timed
    theta = crit(4)
    epsi = r8prem()*pk
    dmgmi = 1.d0-(1.d0+pkc*timef)**d13
!
    ndim = 3
!
    gh1 = yd(8)
    gh2 = yd(9)
    dh1 = dy(8)
    dh2 = dy(9)
    h1 = gh1+theta*dh1
    h2 = gh2+theta*dh2
    h = h1+h2
    dmg = yd(10)
    ddmg = dy(10)
    d = dmg+theta*ddmg
    dp = dy(7)
    dm1 = 1.d0-d
!
!  INITIALISATION DE LA MATRICE DRDY
    call r8inir(nr*nr, 0.d0, drdy, 1)
    do i = 1, 10
        drdy(i, i) = 1.d0
    end do
!
!
!------------CALCUL DES INVARIANTS DE CONTRAINTE  -------
!     attention FGEQUI ne prend pas en compte les SQRT(2)
    do itens = 1, ndt
        epsef(itens) = yd(itens)+theta*dy(itens)
    end do
    call lcopli('ISOTROPE', mod, coefel, hookf)
    sigf(1:ndt) = matmul(hookf(1:ndt, 1:ndt), epsef(1:ndt))
    sigf(1:ndt) = dm1*sigf(1:ndt)
!
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
    seq = equi(1)
    seq0 = seq/(1.d0-d)
    trsig0 = trsig/(1.d0-d)
    if (sequid .eq. 0.d0) then
        sequi = grj0
    else
        sequi = trsig
    end if
!------------ CALCUL DU TENSEUR DEVIATORIQUE DES CONTRAINTES ---
    do itens = 1, ndt
        if (itens .le. 3) then
            dev(itens) = sigf(itens)-trsig/3.d0
        else
            dev(itens) = sigf(itens)*sqrt(2.0d0)
        end if
    end do
!
    shmax = 50.d0
!
    terme1 = (seq*(1.d0-h))/(pk*(1.d0-dmgmi)*(1.d0-d))
    if (seq .le. epsi) then
        sinh1 = 0.d0
    else if (abs(terme1) .lt. shmax) then
        sinh1 = sinh(terme1)
    else
        iret = 1
        goto 9999
    end if
    cosh1 = sqrt(1.d0+sinh1*sinh1)
!
!----- EQUATION DONNANT LA DERIVEE DE L ENDOMMAGEMENT
!
    if (sequi .ge. 0.d0) then
        sinn = alphad*sequi+((1.d0-alphad)*seq)
    else
        sinn = (1.d0-alphad)*seq
    end if
    if ((sinn/sig0) .lt. shmax) then
        sinh2 = sinh(sinn/sig0)
    else
        iret = 1
        goto 9999
    end if
    cosh2 = sqrt(1.d0+sinh2*sinh2)
!
    if (seq .gt. 0.d0) then
!        dFe/dEel
        n(1:ndt) = (1.5d0/seq)*dev(1:ndt)
        call lcprte(n, n, nxn)
        dfedee(1:ndt, 1:ndt) = 1.5d0*id(1:ndt, 1:ndt)
        dfedee(1:ndt, 1:ndt) = dfedee(1:ndt, 1:ndt)-nxn(1:ndt, 1:ndt)
        do i = 1, 6
            do j = 1, 6
                drdy(i, j) = drdy(i, j)+deuxmu*theta*dp/seq*dm1*dfedee(i, j)
            end do
        end do
!        dFe/dp
        do i = 1, 6
            drdy(i, 7) = n(i)
        end do
!        dSeq/dEel
        do i = 1, 6
            dseqde(i) = deuxmu*theta*(1.d0-d)*n(i)
        end do
!        dFp/dEel=
        coef = -eps0*dt*cosh1*(1.d0-h)/pk/(1.d0-dmgmi)/(1.d0-d)
        do i = 1, 6
            drdy(7, i) = coef*dseqde(i)
        end do
!        dFp/dp=
        drdy(7, 7) = 1.d0
!        dFp/dH1=
        drdy(7, 8) = eps0*dt*cosh1*theta*seq0/pk/(1.d0-dmgmi)
        drdy(7, 9) = drdy(7, 8)
!        dFp/dD=
        drdy(7, 10) = 0.d0
!
!        dFH1/dEel=
        coef = ph1*dp*(h1st-delta1*h1)/seq/seq
        do i = 1, 6
            drdy(8, i) = coef*dseqde(i)
        end do
!        dFH2/ds=
        coef = ph2*dp*(h2st-delta2*h2)/seq/seq
        do i = 1, 6
            drdy(9, i) = coef*dseqde(i)
        end do
!        dFH1/dp=
        drdy(8, 7) = -ph1*(h1st-delta1*h1)/seq
        drdy(9, 7) = -ph2*(h2st-delta2*h2)/seq
!        dFH1/dH1=
        drdy(8, 8) = 1.d0+ph1*dp*delta1*theta/seq
        drdy(9, 9) = 1.d0+ph2*dp*delta2*theta/seq
!        dFH1/dD= - ou + ?
        drdy(8, 10) = -ph1*dp*(h1st-delta1*h1)*seq0*theta/seq/seq
        drdy(9, 10) = -ph2*dp*(h2st-delta2*h2)*seq0*theta/seq/seq
!
!        dFD/dEe
        coef = -biga*dt*cosh2/sig0
        do i = 1, 6
            drdy(10, i) = coef*dseqde(i)*(1.d0-alphad)
        end do
!
!        DFD/DD
        drdy(10, 10) = 1.d0+biga*dt*cosh2/sig0*(1-alphad)*theta*seq0
        if (sequid .gt. 0) then
            if (trsig .gt. epsi) then
                drdy(10, 10) = drdy(10, 10)+biga*dt*cosh2/sig0*alphad*theta*trsig0
                do i = 1, 6
                    drdy(10, i) = drdy(10, i)+coef*troisk*theta*alphad*kron(i)*(1.d0-d)
                end do
            end if
        end if
!
    end if
!
!
9999 continue
!
!
end subroutine
