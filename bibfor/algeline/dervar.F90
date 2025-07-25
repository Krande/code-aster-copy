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

subroutine dervar(gamp, nbmat, mater, parame, derpar)
!
    implicit none
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: nbmat
    real(kind=8) :: mater(nbmat, 2), gamp, parame(5), derpar(4)
! --- BUT : CALCUL DES DERIVEES DES VARIABLES D'ECROUISSAGE ------------
! ======================================================================
! IN  : GAMP   : DEFORMATION DEVIATOIRE PLASTIQUE CUMULEE --------------
! --- : NBMAT  : NOMBRE DE PARAMETRES DU MODELE ------------------------
! --- : MATER  : PARAMETRES DU MODELE ----------------------------------
! --- : PARAME : VARIABLES D'ECROUISSAGE (S,A,K,M,O) -------------------
! OUT : DERPAR : DERIVE DES VARIABLES D'ECROUISSAGE PAR RAPPORT A GAMP -
! ------------ : (DS/DGAMP,DA/DGAMP,DK/DGAMP,DM/DGAMP,DO/DGAMP) --------
! ======================================================================
    real(kind=8) :: mun, zero, un, deux, trois
    real(kind=8) :: gamult, gammae, me, ae, mpic, apic, eta, sigc, sigp1, sigp2
    real(kind=8) :: agamp, omega
    real(kind=8) :: ds, domega, dado, da, dk, dmds, dmda, dm
    real(kind=8) :: fact1, fact2, fact3, fact4, puis1
! ======================================================================
! --- INITIALISATION DE PARAMETRES -------------------------------------
! ======================================================================
    parameter(mun=-1.0d0)
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    call jemarq()
! ======================================================================
! --- RECUPERATION DES PARAMETRES DU MODELE ----------------------------
! ======================================================================
    gamult = mater(1, 2)
    gammae = mater(2, 2)
    me = mater(4, 2)
    ae = mater(5, 2)
    mpic = mater(6, 2)
    apic = mater(7, 2)
    eta = mater(8, 2)
    sigc = mater(9, 2)
    sigp1 = mater(13, 2)
    sigp2 = mater(14, 2)
! ======================================================================
! --- RECUPERATION DES PARAMETRES DU MODELE ----------------------------
! ======================================================================
    agamp = parame(2)
    omega = parame(5)
! ======================================================================
! --- CALCUL DES DERIVEES DES VARIABLES D'ECROUISSAGES -----------------
! --- POUR LE CAS GAMP > GAMULT(1-EPS) ---------------------------------
! ======================================================================
    if (gamp .gt. gamult) then
        ds = zero
        da = zero
        dk = zero
        dm = zero
! ======================================================================
! SINON ----------------------------------------------------------------
! ======================================================================
    else
! ======================================================================
! --- CALCUL DE DS/DGAMP = -1/GAMMAE -----------------------------------
! ======================================================================
        if (gamp .lt. gammae) then
            ds = mun/gammae
        else
            ds = zero
        end if
! ======================================================================
! --- CALCUL DE DOMEGA/DGAMP = -----------------------------------------
! ------- (GAMULT-GAMMAE)/(GAMMAE)**ETA*((AE-APIC)/(1-AE))* ------------
! ---- (ETA*GAMP**(ETA-1)/(GAMULT-GAMP) + GAMP**ETA/(GAMULT-GAMP)**2 ---
! ======================================================================
        fact1 = (gamult-gammae)/(gammae**eta)
        fact2 = (ae-apic)/(un-ae)
        fact3 = eta*(gamp**(eta-un))*(gamult-gamp)
        fact4 = gamp**eta
        domega = fact1*fact2*(fact3+fact4)
! ======================================================================
! --- CALCUL DE DA/DOMEGA = (1-APIC)/(1+OMEGA(GAMP))**2 ----------------
! ======================================================================
        dado = (un-apic)/((gamult-gamp+omega)*(gamult-gamp+omega))
! ======================================================================
! --- CALCUL DE DA/DGAMP = DA/DOMEGA * DOMEGA/DGAMP --------------------
! ======================================================================
        da = dado*domega
! ======================================================================
! --- CALCUL DE DK/DGAMP = -(2/3)**(1/(2*A(GAMP)))* --------------------
! ------------------------  LOG(2/3)/(2*A(GAMP)*A(GAMP))*DA/DGAMP ------
! ======================================================================
        puis1 = un/(deux*agamp)
        dk = mun*((deux/trois)**puis1)*log(deux/trois)*da/(deux*agamp*agamp)
! ======================================================================
! --- CALCUL DE DM/DGAMP = ---------------------------------------------
! --- SI GAMP  <  GAMMAE : DM/DA*DA/DGAMP+DM/DS*DS/DGAMP ---------------
! --- AVEC DM/DA = -SIGC/SIGP1*LOG(MPIC*SIGP1/SIGC+1)*APIC/A(GAMP)**2* -
! ---------------- (MPIC*SIGP1/SIGC+1)**(APIC/A(GAMP)) -----------------
! -------- DM/DS = -SIGC/SIGP1 -----------------------------------------
! --- SI GAMP  >= GAMMAE : DM/DA*DA/DGAMP ------------------------------
! --- AVEC DM/DA = -SIGC/SIGP2*LOG(ME*SIGP2/SIGC)*AE/A(GAMP)**2* -------
! ---------------- (ME*SIGP2/SIGC)**(AE/A(GAMP)) -----------------------
! ======================================================================
        if (gamp .lt. gammae) then
            dmds = mun*sigc/sigp1
            fact1 = mpic*sigp1/sigc+un
            fact2 = apic/agamp
            dmda = mun*sigc*log(fact1)*fact2*(fact1**fact2)/(sigp1*agamp)
            dm = dmda*da+dmds*ds
        else
            fact1 = me*sigp2/sigc
            fact2 = ae/agamp
            dmda = mun*sigc*log(fact1)*fact2*(fact1**fact2)/(sigp2*agamp)
            dm = dmda*da
        end if
    end if
! ======================================================================
! --- STOCKAGE ---------------------------------------------------------
! ======================================================================
    derpar(1) = ds
    derpar(2) = da
    derpar(3) = dk
    derpar(4) = dm
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
