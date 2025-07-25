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
subroutine dpvpre(mod, nvi, option, crit, instam, &
                  instap, nbmat, materf, sigm, deps, &
                  vim, vip, sig, nbre, dsidep, &
                  iret)
! =====================================================================
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dpvpcr.h"
#include "asterfort/dpvpdb.h"
#include "asterfort/dpvpot.h"
#include "asterfort/dpvpsi.h"
#include "asterfort/dpvpva.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
#include "asterfort/trace.h"
    integer(kind=8) :: iret, nvi, nbmat
    real(kind=8) :: deps(6), vim(nvi), vip(nvi), sig(6)
    real(kind=8) :: sigm(6), materf(nbmat, 2), dsidep(6, 6), crit(3)
    real(kind=8) :: instam, instap
    character(len=8) :: mod
    character(len=16) :: option
! =====================================================================
! --- LOI DE COMPORTEMENT DRUCKER PRAGER VISCOPLASTIQUE ---------------
! --- VISC_DRUC_PRAG --------------------------------------------------
! ----RESOLUTION -----------------------------------------------------
    aster_logical :: resi
    integer(kind=8) :: ndt, ndi, ii, pos, nbre
    real(kind=8) :: deux, trois
    real(kind=8) :: ppic, pult
    real(kind=8) :: dp
    real(kind=8) :: dt, plas
    real(kind=8) :: fcrit
    real(kind=8) :: siim, siie, seqm, seqe, i1m, i1e
    real(kind=8) :: fonecm(3), fonecp(3)
    real(kind=8) :: hookf(6, 6), dkooh(6, 6)
    real(kind=8) :: sige(6), se(6), sm(6)
    real(kind=8) :: inte(6)
! =====================================================================
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
    common/tdim/ndt, ndi
! =====================================================================
! --- RECUPERATION DES MATERIAUX --------------------------------------
! =====================================================================
    ppic = materf(4, 2)
    pult = materf(5, 2)
! =====================================================================
! --- INITIALISATION --------------------------------------------------
! =====================================================================
    hookf(:, :) = 0.0d0
    dkooh(:, :) = 0.0d0
    sig(:) = 0.0d0
    sige(:) = 0.0d0
!
    iret = 0
    dt = instap-instam
    plas = 0.d0
    dp = 0.d0
!
! =====================================================================
! =====================================================================
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
    ASSERT((option(1:9) .eq. 'RIGI_MECA') .or. resi)
! =====================================================================
! --- OPERATEUR ELASTIQUE LINEAIRE ISOTROPE ---------------------------
! =====================================================================
    call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
    call lcopil('ISOTROPE', mod, materf(1, 1), dkooh)
!
    call lcdevi(sigm, sm)
    siim = dot_product(sm(1:ndt), sm(1:ndt))
    seqm = sqrt(trois*siim/deux)
    i1m = trace(ndi, sigm)
!
! =====================================================================
! --- PREDICTION ELASTIQUE : SIGF = HOOKF EPSP -----------------------
! =====================================================================
!
!
    inte(1:ndt) = matmul(hookf(1:ndt, 1:ndt), deps(1:ndt))
    sige(1:ndt) = sigm(1:ndt)+inte(1:ndt)
!
    call lcdevi(sige, se)
    siie = dot_product(se(1:ndt), se(1:ndt))
    seqe = sqrt(trois*siie/deux)
    i1e = trace(ndi, sige)
! =====================================================================
! --- CALCUL DU CRITERE -----------------------------------------------
! =====================================================================
    if (resi) then
! =====================================================================
! --- SIGNE DU CRITERE ------------------------------------------------
! =====================================================================
        call dpvpva(vim(1), nbmat, materf, fonecm)
!
        fcrit = dpvpcr(fonecm, seqe, i1e)
!
        if (fcrit .gt. 0.0d0) then
            plas = 1.0d0
! =====================================================================
! ---------------------RESOLUTION EQ NON LINEAIRE EN DP----------------
! =====================================================================
            call dpvpdb(nbmat, materf, crit, dt, vim, &
                        vip, nvi, seqe, i1e, seqm, &
                        i1m, dp, nbre, iret)
        else
            plas = 0.0d0
            dp = 0.0d0
            nbre = 0
        end if
! =====================================================================
! --- MISE A JOUR DES CONTRAINTES TENANT COMPTE DE DP SI VISCOPLASTICIT
! =====================================================================
        if (plas .eq. 0.0d0) then
            do ii = 1, ndt
                sig(ii) = sige(ii)
            end do
!
            vip(1) = vim(1)
            vip(3) = vim(3)
            vip(4) = vim(4)
        else
            vip(1) = vim(1)+dp
!
            call dpvpva(vip(1), nbmat, materf, fonecp)
            call dpvpsi(nbmat, materf, se, seqe, i1e, &
                        fonecp, dp, sig)
        end if
! =====================================================================
! --- STOCKAGE DES VARIABLES INTERNES ---------------------------------
! =====================================================================
        vip(2) = plas
!
        vip(nvi) = nbre
!
        if (vip(1) .lt. ppic) then
            pos = 1
        else if (vip(1) .lt. pult) then
            pos = 2
        else
            pos = 3
        end if
!
        vip(3) = pos
    end if
! =====================================================================
! --- TERMES DE L OPERATEUR TANGENT -----------------------------------
! =====================================================================
    if (option(10:14) .eq. '_ELAS') then
        dsidep(1:ndt, 1:ndt) = hookf(1:ndt, 1:ndt)
    end if
!
    if (option(1:14) .eq. 'RIGI_MECA_TANG') then
        dsidep(1:ndt, 1:ndt) = hookf(1:ndt, 1:ndt)
    end if
    if (option(1:9) .eq. 'FULL_MECA') then
        if (vip(2) .eq. 1.d0) then
            call dpvpot(mod, vim(1), vip(1), nbmat, materf, &
                        sige, dt, dp, plas, dsidep)
        else if (vip(2) .eq. 0.d0) then
            dsidep(1:ndt, 1:ndt) = hookf(1:ndt, 1:ndt)
        end if
    end if
! =====================================================================
end subroutine
