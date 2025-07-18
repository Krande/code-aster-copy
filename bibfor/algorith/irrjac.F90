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

subroutine irrjac(fami, kpg, ksp, mod, nmat, &
                  mater, yf, dy, nmod, drdy)
!
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/irrfss.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcicma.h"
#include "asterfort/lcnrts.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
#include "asterfort/rcvarc.h"
    character(len=*) :: fami
    character(len=8) :: mod
    integer(kind=8) :: nmat, nmod, kpg, ksp
    real(kind=8) :: mater(nmat, 2), yf(*), dy(*), drdy(nmod, nmod)
! person_in_charge: jean-luc.flejou at edf.fr
!
!     ------------------------------------------------------------------
!  IRRAD3M    : CALCUL DU JACOBIEN DU SYSTEME NL A RESOUDRE = DRDY
!               DY    = ( DSIG  DP  DETA   DPI   DG  (DEPS3) )
!               Y     = ( SIG   P   ETA    PI    G   (DEPS3) )
!               DRDY  = ( DRSDS  DRSDP  DRSDE  DRSDI DRSDG (DRSDE3) )
!                       ( DRPDS  DRPDP  DRPDE  DRPDI DRPDG (DRPDE3) )
!                       ( DREDS  DREDP  DREDE  DREDI DREDG (DREDE3) )
!                       ( DRIDS  DRIDP  DRIDE  DRIDI DRIDG (DRIDE3) )
!                       ( DRGDS  DRGDP  DRGDE  DRGDI DRGDG (DRGDE3) )
!                       ((DQDS) (DQDP) (DQDE) (DQDI) (DQDG) (DQDE3) )
!  IN  FAMI   :  FAMILLE DES POINTS DE GAUSS
!      KPG    :  NUMERO DU POINT DE GAUSS
!      KSP    :  NUMERO DU SOUS POINT DE GAUSS
!      MOD    :  TYPE DE MODELISATION
!      NMAT   :  DIMENSION MATER
!      MATER  :  COEFFICIENTS MATERIAU A T+DT
!      YF     :  VARIABLES A T + DT
!                ( SIGF  PF   ETAF   PIF   GF  (EPS3F) )
!      DY     :  SOLUTION
!                ( DSIG  DP   DETA   DPI   SG  (DEPS3) )
!      NMOD   :  DIMENSION DECLAREE DRDY
!  OUT DRDY   :  JACOBIEN DU SYSTEME NON LINEAIRE
!  ----------------------------------------------------------------
    real(kind=8) :: fkooh(6, 6), k, n, p0, ai0, irrad, irraf, sigf(6)
    real(kind=8) :: etais, pf, dp, dpi, dphi, dev(6), dfds(6), zetaf
    real(kind=8) :: detai, seqf
    real(kind=8) :: drsds(6, 6), drsdp(6), drsde(6), drsdi(6), drsdg(6)
    real(kind=8) :: drpds(6), drpdp, drpde, drpdi, drpdg
    real(kind=8) :: dreds(6), dredp, drede, dredi, dredg
    real(kind=8) :: drids(6), dridp, dride, dridi, dridg
    real(kind=8) :: drgds(6), drgdp, drgde, drgdi, drgdg
    real(kind=8) :: drsde3(6), drpde3, drede3, dride3, drgde3
    real(kind=8) :: dqds(4), dqdp, dqde, dqdi, dqdg, dqde3
    real(kind=8) :: sr, ddfdds(6, 6), etaif, dede3(6), hookf(6, 6)
    real(kind=8) :: pk, kappa, r02, pe, penpe, spe
    integer(kind=8) :: ndt, ndi, iret
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
    data dede3/0.0d0, 0.0d0, -1.0d0, 0.0d0, 0.0d0, 0.0d0/
!
    call lcopil('ISOTROPE', mod, mater(1, 1), fkooh)
    call lcopli('ISOTROPE', mod, mater(1, 1), hookf)
    call rcvarc('F', 'IRRA', '-', fami, kpg, &
                ksp, irrad, iret)
    call rcvarc('F', 'IRRA', '+', fami, kpg, &
                ksp, irraf, iret)
!     ARRET DANS IRRMAT SI  IRRAD .GT. IRRAF*1.00001
    if (irrad .gt. irraf) then
        dphi = 0.0d0
    else
        dphi = irraf-irrad
    end if
!
!     RECUPERATION DES INCREMENTS DES VARIABLES INTERNES
    dp = dy(ndt+1)
    detai = dy(ndt+2)
    dpi = dy(ndt+3)
!
!     RECUPERATION DES VARIABLES INTERNES A t+
    sigf(1:ndt) = yf(1:ndt)
    pf = yf(ndt+1)
    etaif = yf(ndt+2)
!
!     CARACTERISTIQUES MATERIAUX
    ai0 = mater(4, 2)
    etais = mater(5, 2)
    k = mater(7, 2)
    n = mater(8, 2)
    p0 = mater(9, 2)
    kappa = mater(10, 2)
    r02 = mater(11, 2)
    zetaf = mater(12, 2)
    penpe = mater(13, 2)
    pk = mater(14, 2)
    pe = mater(15, 2)
    spe = mater(16, 2)
!
    call lcdevi(sigf, dev)
    seqf = lcnrts(dev)
    if (seqf .eq. 0.0d0) then
        dfds(:) = 0.0d0
    else
        dfds(1:ndt) = (1.5d0/seqf)*dev(1:ndt)
    end if
!
! - DRSDS
    call irrfss(sigf, ddfdds)
    drsds(1:ndt, 1:ndt) = (dp+dpi)*ddfdds(1:ndt, 1:ndt)
    drsds(1:ndt, 1:ndt) = fkooh(1:ndt, 1:ndt)+drsds(1:ndt, 1:ndt)
! - DRSDP
    drsdp(1:ndt) = dfds(1:ndt)
! - DRSDE
    drsde(:) = 0.0d0
! - DRSDI
    drsdi(1:ndt) = dfds(1:ndt)
! - DRSDG
!       CALL LCEQVN(NDT,ID,DRSDG)
    drsdg(:) = 0.0d0
!
! - LOI DE COMPORTEMENT
    if (pf .lt. pk) then
        sr = kappa*r02
    else if (pf .lt. pe) then
        sr = penpe*(pf-pe)+spe
    else
        sr = k*((pf+p0)**n)
    end if
! - DRPDS
    if (((seqf .ge. sr) .and. (dp .ge. 0.0d0)) .or. (dp .gt. r8prem())) then
        drpds(1:ndt) = dfds(1:ndt)
        drpds(1:ndt) = (1.0d0/hookf(1, 1))*drpds(1:ndt)
    else
        drpds(:) = 0.0d0
    end if
! - DRPDP
    if (((seqf .ge. sr) .and. (dp .ge. 0.0d0)) .or. (dp .gt. r8prem())) then
        if (pf .lt. pk) then
            drpdp = 0.0d0
        else if (pf .lt. pe) then
            drpdp = -penpe/hookf(1, 1)
        else
            drpdp = (-n*k*((pf+p0)**(n-1.0d0)))/hookf(1, 1)
        end if
    else
        drpdp = 1.0d0
    end if
! - DRPDE
    drpde = 0.0d0
! - DRPDI
    drpdi = 0.0d0
! - DRPDG
    drpdg = 0.0d0
!
! - DREDS
    dreds(1:ndt) = (-dphi*zetaf*0.50d0/hookf(1, 1))*dfds(1:ndt)
! - DREDP
    dredp = 0.0d0
! - DREDE
    drede = 1.0d0/hookf(1, 1)
! - DREDI
    dredi = 0.0d0
! - DREDG
    dredg = 0.0d0
!
! - DRIDS
    drids(:) = 0.0d0
! - DRIDP
    dridp = 0.0d0
! - DRIDE
    if ((etaif-detai) .gt. etais) then
        dride = -ai0
    else if (etaif .le. etais) then
        dride = 0.0d0
    else
        dride = -ai0
    end if
!
! - DRIDI
    dridi = 1.0d0
! - DRIDG
    dridg = 0.0d0
!
!
! - DRGDS
    drgds(:) = 0.0d0
! - DRGDP
    drgdp = 0.0d0
! - DRGDE
    drgde = 0.0d0
! - DRGDI
    drgdi = 0.0d0
! - DRGDG
    drgdg = 1.0d0
!
! - CONTRAINTES PLANES
    if (mod(1:6) .eq. 'C_PLAN') then
! - DRSDE3
        drsde3(1:ndt) = dede3(1:ndt)
! - DRPDE3
        drpde3 = 0.0d0
! - DREDE3
        drede3 = 0.0d0
! - DRIDE3
        dride3 = 0.0d0
! - DRGDE3
        drgde3 = 0.0d0
! - DQDS
        dqds(1) = (-(dp+dpi)*(hookf(3, 3)*ddfdds(3, 1)+hookf(3, 1)* &
                              ddfdds(1, 1)+hookf(3, 2)*ddfdds(2, 1)+hookf(3, 4)*ddfdds(4, 1)) &
                   )/hookf(1, 1)
        dqds(2) = (-(dp+dpi)*(hookf(3, 3)*ddfdds(3, 2)+hookf(3, 1)* &
                              ddfdds(1, 2)+hookf(3, 2)*ddfdds(2, 2)+hookf(3, 4)*ddfdds(4, 2)) &
                   )/hookf(1, 1)
        dqds(3) = (-(dp+dpi)*(hookf(3, 3)*ddfdds(3, 3)+hookf(3, 1)* &
                              ddfdds(1, 3)+hookf(3, 2)*ddfdds(2, 3)+hookf(3, 4)*ddfdds(4, 3)) &
                   )/hookf(1, 1)
        dqds(4) = (-(dp+dpi)*(hookf(3, 3)*ddfdds(3, 4)+hookf(3, 1)* &
                              ddfdds(1, 4)+hookf(3, 2)*ddfdds(2, 4)+hookf(3, 4)*ddfdds(4, 4)) &
                   )/hookf(1, 1)
! - DQDP
        dqdp = ( &
               -hookf(3, 1)*dfds(1)-hookf(3, 2)*dfds(2)-hookf(3, 3)*dfds(3)-hookf(3, 4)*d&
               &fds(4))/hookf(1, &
               1 &
               )
! - DQDE
        dqde = 0.0d0
! - DQDI
        dqdi = ( &
               -hookf(3, 1)*dfds(1)-hookf(3, 2)*dfds(2)-hookf(3, 3)*dfds(3)-hookf(3, 4)*d&
               &fds(4))/hookf(1, &
               1 &
               )
! - DQDG
!        DQDG=-HOOKF(3,3)/HOOKF(1,1)
        dqdg = 0.0d0
! - DQDE3
        dqde3 = hookf(3, 3)/hookf(1, 1)
    end if
!
! - ASSEMBLAGE
!
! - DRDY(T+DT)  =  DRSDS  DRSDP  DRSDE  DRSDI DRSDG (DRSDE3)
!                  DRPDS  DRPDP  DRPDE  DRPDI DRPDG (DRPDE3)
!                  DREDS  DREDP  DREDE  DREDI DREDG (DREDE3)
!                  DRIDS  DRIDP  DRIDE  DRIDI DRIDG (DRIDE3)
!                  DRGDS  DRGDP  DRGDE  DRGDI DRGDG (DRGDE3)
!                 (DQDS) (DQDP) (DQDE) (DQDI) (DQDG) (DQDE3)
!
!
    call lcicma(drsds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                1, 1)
    call lcicma(drsdp, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+1)
    call lcicma(drsde, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+2)
    call lcicma(drsdi, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+3)
    call lcicma(drsdg, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+4)
!
    call lcicma(drpds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+1, 1)
    call lcicma(dreds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+2, 1)
    call lcicma(drids, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+3, 1)
    call lcicma(drgds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+4, 1)
!
    drdy(ndt+1, ndt+1) = drpdp
    drdy(ndt+1, ndt+2) = drpde
    drdy(ndt+1, ndt+3) = drpdi
    drdy(ndt+1, ndt+4) = drpdg
!
    drdy(ndt+2, ndt+1) = dredp
    drdy(ndt+2, ndt+2) = drede
    drdy(ndt+2, ndt+3) = dredi
    drdy(ndt+2, ndt+4) = dredg
!
    drdy(ndt+3, ndt+1) = dridp
    drdy(ndt+3, ndt+2) = dride
    drdy(ndt+3, ndt+3) = dridi
    drdy(ndt+3, ndt+4) = dridg
!
    drdy(ndt+4, ndt+1) = drgdp
    drdy(ndt+4, ndt+2) = drgde
    drdy(ndt+4, ndt+3) = drgdi
    drdy(ndt+4, ndt+4) = drgdg
!
    if (mod(1:6) .eq. 'C_PLAN') then
!
        call lcicma(drsde3, 6, 1, ndt, 1, &
                    1, 1, drdy, nmod, nmod, &
                    1, ndt+5)
        call lcicma(dqds, 1, 6, 1, ndt, &
                    1, 1, drdy, nmod, nmod, &
                    ndt+5, 1)
        drdy(ndt+1, ndt+5) = drpde3
        drdy(ndt+2, ndt+5) = drede3
        drdy(ndt+3, ndt+5) = dride3
        drdy(ndt+4, ndt+5) = drgde3
        drdy(ndt+5, ndt+1) = dqdp
        drdy(ndt+5, ndt+2) = dqde
        drdy(ndt+5, ndt+3) = dqdi
        drdy(ndt+5, ndt+4) = dqdg
        drdy(ndt+5, ndt+5) = dqde3
!
    end if
end subroutine
