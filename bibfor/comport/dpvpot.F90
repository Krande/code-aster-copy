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
subroutine dpvpot(mod, vim, vip, nbmat, mater, &
                  sig, dt, dp, plas, dsidep)
!
    implicit none
#include "asterfort/dpvpdv.h"
#include "asterfort/dpvpva.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcprte.h"
#include "asterfort/trace.h"
    integer(kind=8) :: ndt, ndi
    integer(kind=8) :: nbmat
    real(kind=8) :: dt, dp, plas
    real(kind=8) :: mater(nbmat, 2), vim(4), vip(4), sig(6), dsidep(6, 6)
    character(len=8) :: mod
! --- BUT   OPERATEUR TANGENT COHERENT POUR LA LOI --------------------
! --- VISC_DRUC_PRAG --------------------------------------------------
! =====================================================================
    integer(kind=8) :: ii, jj
    real(kind=8) :: zero, un, deux, trois, neuf, unstr
    real(kind=8) :: troisk, deuxmu, k, mu
    real(kind=8) :: pref, a, n
    real(kind=8) :: fonc1, fonc2, fonc3, fonc4, fonc, foncp
    real(kind=8) :: beta
    real(kind=8) :: alpham, betam, rm
    real(kind=8) :: dalpdp, dbetdp, drdp
    real(kind=8) :: fonecm(3), fonecp(3), fonder(3)
    real(kind=8) :: sii, seq, i1
    real(kind=8) :: scal1, scal3, scal4, scal5, scal6
    real(kind=8) :: scal12
    real(kind=8) :: denom, dfdp, const, const1
    real(kind=8) :: kron(6)
    real(kind=8) :: dsede(6, 6)
    real(kind=8) :: s(6)
    real(kind=8) :: dsdsig(6, 6), dqdsig(6), dfdsig(6), dpdsig(6)
    real(kind=8) :: adidsi(6), bdidsi(6), cdidsi(6)
    real(kind=8) :: dgdsig(6)
    real(kind=8) :: dqdeps(6, 6), dsdeps(6, 6), dpdeps(6)
    real(kind=8) :: di1ede(6), di1de(6)
    real(kind=8) :: vect1(6), vect2(6), vect3(6), vect4(6)
    real(kind=8) :: matr1(6, 6), matr2(6, 6), matr3(6, 6)
    real(kind=8) :: matr1a(6, 6), matr1b(6, 6)
    real(kind=8) :: part1(6, 6), part2(6, 6), part3(6, 6), part4(6, 6)
    real(kind=8) :: inter1(6, 6), inter2(6, 6), inter3(6, 6)
    real(kind=8) :: int2a(6, 6), int2b(6, 6), dsdept(6, 6)
    real(kind=8) :: tol
! =====================================================================
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
    parameter(neuf=9.0d0)
    parameter(unstr=1.d0/3.0d0)
    parameter(tol=1.d-12)
!
! =====================================================================
    common/tdim/ndt, ndi
! =================================================================
    data kron/un, un, un, zero, zero, zero/
! =================================================================
! ---- RECUPERATION DES PARAMETRES MATERIAUX ----------------------
! =================================================================
!
    mu = mater(4, 1)
    k = mater(5, 1)
    pref = mater(1, 2)
    a = mater(2, 2)
    n = mater(3, 2)
    troisk = trois*k
    deuxmu = deux*mu
! =====================================================================
! --- INITIALISATIONS DES VECTEURS ------------------------------------
! =====================================================================
    vect1(:) = 0.0d0
    vect2(:) = 0.0d0
    vect3(:) = 0.0d0
    vect4(:) = 0.0d0
    dqdsig(:) = 0.0d0
    dfdsig(:) = 0.0d0
    dpdsig(:) = 0.0d0
    dgdsig(:) = 0.0d0
    dpdeps(:) = 0.0d0
    di1ede(:) = 0.0d0
    di1de(:) = 0.0d0
    adidsi(:) = 0.0d0
    bdidsi(:) = 0.0d0
    cdidsi(:) = 0.0d0
    s(:) = 0.0d0
! =====================================================================
! --- INITIALISATIONS DES MATRICES ------------------------------------
! =====================================================================
    matr1(:, :) = 0.0d0
    matr1a(:, :) = 0.0d0
    matr1b(:, :) = 0.0d0
    matr2(:, :) = 0.0d0
    matr3(:, :) = 0.0d0
    part1(:, :) = 0.0d0
    part2(:, :) = 0.0d0
    part3(:, :) = 0.0d0
    part4(:, :) = 0.0d0
    inter1(:, :) = 0.0d0
    inter2(:, :) = 0.0d0
    int2a(:, :) = 0.0d0
    int2b(:, :) = 0.0d0
    inter3(:, :) = 0.0d0
    dsdsig(:, :) = 0.0d0
    dsdeps(:, :) = 0.0d0
    dsdept(:, :) = 0.0d0
    dqdeps(:, :) = 0.0d0
    dsidep(:, :) = 0.0d0
    dsede(:, :) = 0.0d0
!
    call lcopli('ISOTROPE', mod, mater(1, 1), dsede)
! =====================================================================
! --- CAS ELASTIQUE ---------------------------------------------------
    if ((plas .eq. 0.0d0) .or. (dp .eq. 0.d0) .or. (abs(dp) .lt. tol)) then
        dsidep(1:ndt, 1:ndt) = dsede(1:ndt, 1:ndt)
        goto 999
    else
! =================================================================
! ----  CALCUL DU DEVIATEUR - DE LA CONTRAINTE EQUIVALENTE  -------
! ----  ET DE LA TRACE --------------------------------------------
! =================================================================
        call lcdevi(sig, s)
        sii = dot_product(s(1:ndt), s(1:ndt))
        seq = sqrt(trois*sii/deux)
        i1 = trace(ndi, sig)
!
! =====================================================================
! --- FONCTIONS D ECROUISSAGE ET LEURS DERIVEES------------------------
! =====================================================================
        call dpvpva(vim, nbmat, mater, fonecm)
        call dpvpva(vip, nbmat, mater, fonecp)
        call dpvpdv(vip, nbmat, mater, fonder)
!
        alpham = fonecm(1)
        rm = fonecm(2)
        betam = fonecm(3)
!
        beta = fonecp(3)
!
        dalpdp = fonder(1)
        drdp = fonder(2)
        dbetdp = fonder(3)
!
        const = a*dt/(pref)**n
! =====================================================================
! --- CALCUL DE DSIDEP ------------------------------------------------
! =====================================================================
        do ii = 1, ndi
            do jj = 1, ndi
                dsdsig(ii, jj) = -un/trois
            end do
        end do
        do ii = 1, ndt
            dsdsig(ii, ii) = dsdsig(ii, ii)+un
        end do
!
!
! =====================================================================
! --- CALCUL DE LA TROISIEME PARTIE DU TERME DS/DEPS ------------------
! =====================================================================
! --- CALCUL DE DSIEQ/ DSIG -------------------------------------------
! =====================================================================
!
        scal1 = trois/deux/seq
!
        dqdsig(1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), s(1:ndt))
        dqdsig(1:ndt) = scal1*dqdsig(1:ndt)
! =====================================================================
! --- CALCUL DE ALPHA * DI1/ DSIG -------------------------------------
! =====================================================================
        adidsi(1:ndt) = alpham*kron(1:ndt)
! =====================================================================
! --- CALCUL DE ALPHA_CONS * DI1/ DSIG * DP ---------------------------
! =====================================================================
        scal12 = dalpdp*dp
        bdidsi(1:ndt) = scal12*kron(1:ndt)
! =====================================================================
! --- CALCUL DE Df/ DSIG ----------------------------------------------
! =====================================================================
        cdidsi(1:ndt) = dqdsig(1:ndt)+adidsi(1:ndt)
        dfdsig(1:ndt) = cdidsi(1:ndt)+bdidsi(1:ndt)
! =====================================================================
! --- CALCUL DE DfDp --------------------------------------------------
! =====================================================================
        fonc1 = seq+alpham*i1-rm
!
        fonc2 = trois*mu+drdp-dalpdp*i1+neuf*k*alpham*betam
!
        fonc3 = neuf*k*(alpham*dbetdp+betam*dalpdp)
!
        fonc4 = neuf*k*dalpdp*dbetdp
!
!
        fonc = fonc1-fonc2*dp-fonc3*dp**2-fonc4*dp**3
        foncp = -fonc2-deux*dp*fonc3-trois*dp**2*fonc4
!
        if (fonc .gt. zero) then
            fonc = fonc
        else
            dsidep(1:ndt, 1:ndt) = dsede(1:ndt, 1:ndt)
            goto 999
        end if
!
        const1 = n*const*fonc**(n-un)
        dfdp = const1*foncp-un
!
        if (dfdp .eq. zero) then
            dsidep(1:ndt, 1:ndt) = dsede(1:ndt, 1:ndt)
            goto 999
        else
            denom = -un/dfdp
        end if
! =====================================================================
! --- CALCUL DE Df/ DSIG ----------------------------------------------
! =====================================================================
        dfdsig(1:ndt) = const1*dfdsig(1:ndt)
! =====================================================================
! --- CALCUL DE d deltap/dSIG -----------------------------------------
! =====================================================================
        dpdsig(1:ndt) = denom*dfdsig(1:ndt)
! =====================================================================
! --- CALCUL DE d deltap/dEPS -----------------------------------------
! =====================================================================
        dpdeps(1:ndt) = matmul(dsede(1:ndt, 1:ndt), dpdsig(1:ndt))
! =====================================================================
! --- CALCUL DE 3GDT/SEQ *se * deltap/dEPS ----------------------------
! =====================================================================
        call lcprte(dpdeps, s, matr1a)
! =====================================================================
! --- TRANSPOSEE ------------------------------------------------------
! =====================================================================
        matr1b(1:ndt, 1:ndt) = transpose(matr1a(1:ndt, 1:ndt))
! =====================================================================
! --- SYMETRISATION  --------------------------------------------------
! =====================================================================
        do ii = 1, ndt
            do jj = 1, ndt
                matr1(ii, jj) = un/deux*(matr1a(ii, jj)+matr1b(ii, jj))
            end do
        end do
!
        scal3 = -trois*mu/seq
!
        part3(1:ndt, 1:ndt) = scal3*matr1(1:ndt, 1:ndt)
!
! =====================================================================
! --- CALCUL DE LA  PREMIERE PARTIE DU TERME DS/DEPS ------------------
! =====================================================================
! --- CALCUL DE dse / deps *(1-3GDP/SEQ) ----------------------------
! =====================================================================
        scal4 = deuxmu*(un-trois*mu*dp/seq)
        part1(1:ndt, 1:ndt) = scal4*dsdsig(1:ndt, 1:ndt)
! =====================================================================
! --- CALCUL DE LA  DEUXIEME PARTIE DU TERME DS/DEPS ------------------
! =====================================================================
! --- CALCUL DE 3GDP/SEQ**2 *(se * dSEQ/dEPS ------------------------
! =====================================================================
        scal5 = neuf*mu*mu*dp/seq/seq/seq
        call lcprte(s, s, matr3)
!
        part2(1:ndt, 1:ndt) = scal5*matr3(1:ndt, 1:ndt)
! =====================================================================
! --- SOMMATION DES PARTIES DE ds/dEPS --------------------------------
! =====================================================================
        inter1(1:ndt, 1:ndt) = part1(1:ndt, 1:ndt)+part2(1:ndt, 1:ndt)
!
        dsdeps(1:ndt, 1:ndt) = inter1(1:ndt, 1:ndt)+part3(1:ndt, 1:ndt)
!
        dsdept(1:ndt, 1:ndt) = transpose(dsdeps(1:ndt, 1:ndt))
        do ii = 1, ndt
            do jj = 1, ndt
                dsdeps(ii, jj) = un/deux*(dsdeps(ii, jj)+dsdept(ii, jj))
            end do
        end do
! =====================================================================
! --- CALCUL  DU TERME DI/DEPS ----------------------------------------
! =====================================================================
! --- CALCUL DE dI1E/dEPS ---------------------------------------------
! =====================================================================
        di1ede(1:ndt) = troisk*kron(1:ndt)
! =====================================================================
! --- CALCUL DE 9KBETAdp/dEPS -----------------------------------------
! =====================================================================
        scal6 = -neuf*k*beta
        vect2(1:ndt) = scal6*dpdeps(1:ndt)
! =====================================================================
! --- CALCUL DE dI1/dEPS ----------------------------------------------
! =====================================================================
        di1de(1:ndt) = di1ede(1:ndt)+vect2(1:ndt)
! =====================================================================
! --- CALCUL DE I * dI/dEPS -------------------------------------------
! =====================================================================
        call lcprte(kron, di1de, int2a)
! =====================================================================
! --- TRANSPOSEE DE I * dI/dEPS ---------------------------------------
! =====================================================================
        int2b(1:ndt, 1:ndt) = transpose(int2a(1:ndt, 1:ndt))
! =====================================================================
! --- SYMETRISATION  --------------------------------------------------
! =====================================================================
        do ii = 1, ndt
            do jj = 1, ndt
                inter2(ii, jj) = un/deux*(int2a(ii, jj)+int2b(ii, jj))
            end do
        end do
        inter2(1:ndt, 1:ndt) = unstr*inter2(1:ndt, 1:ndt)
        dsidep(1:ndt, 1:ndt) = dsdeps(1:ndt, 1:ndt)+inter2(1:ndt, 1:ndt)
    end if
! =====================================================================
999 continue
! =====================================================================
end subroutine
