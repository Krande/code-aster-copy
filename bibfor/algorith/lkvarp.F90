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

subroutine lkvarp(vin, nbmat, mater, paraep)
!
    implicit none
    integer(kind=8) :: nbmat
    real(kind=8) :: vin(7), paraep(3), mater(nbmat, 2)
! --- MODELE LETK : LAIGLE VISCOPLASTIQUE-------------------------
! ================================================================
! --- BUT : CALCUL DES VARIABLES D'ECROUISSAGE -------------------
! ================================================================
! --- : VIN    : TABLEAU DES VARIABLE INTERNES (ICI XIP) ---------
! --- : NBMAT  : NOMBRE DE PARAMETRES DU MODELE ------------------
! --- : MATER  : PARAMETRES DU MODELE ----------------------------
! OUT : PARAEP : VARIABLE D'ECROUISSAGE --------------------------
! ------------ : AXIP, SXIP, MXIP --------------------------------
! ================================================================
    real(kind=8) :: sxip, axip, mxip
    real(kind=8) :: xiult, xie, xipic, m0, mult, me, mpic, a0, ae, apic
    real(kind=8) :: s0, xams, eta, sigc, xip
    real(kind=8) :: sigp1, sigp2, un, zero
    real(kind=8) :: fact1, fact2, fact3, fact4
! ================================================================
! --- INITIALISATION DE PARAMETRES -------------------------------
! ================================================================
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
! ================================================================
! --- RECUPERATION DE PARAMETRES DU MODELE -----------------------
! ================================================================
    sigc = mater(3, 2)
    xams = mater(6, 2)
    eta = mater(7, 2)
    a0 = mater(8, 2)
    ae = mater(9, 2)
    apic = mater(10, 2)
    s0 = mater(11, 2)
    m0 = mater(12, 2)
    me = mater(13, 2)
    mpic = mater(14, 2)
    mult = mater(15, 2)
    xiult = mater(16, 2)
    xie = mater(17, 2)
    xipic = mater(18, 2)
!
    sigp1 = mater(23, 2)
!
    sigp2 = ((mult*(sigc)**(ae-un))/(me**ae))**(un/(ae-un))
!
! ================================================================
! CALCUL DES VARIABLES D'ECROUISSAGES POUR LE CAS 0<XIP<XIPIC-----
! ================================================================
    xip = vin(1)
!
!
    if ((xip .ge. zero) .and. (xip .lt. xipic)) then
        fact1 = log(un+xip/xams/xipic)
        fact2 = (apic-a0)/log(un+un/xams)
        axip = a0+fact1*fact2
!
        fact3 = (mpic-m0)/log(un+un/xams)
        mxip = m0+fact1*fact3
!
        fact4 = (un-s0)/log(un+un/xams)
        sxip = s0+fact1*fact4
! ================================================================
! CALCUL DES VARIABLES D'ECROUISSAGES POUR LE CAS XIPIX< XIP < XIE
! ================================================================
    else if ((xip .ge. xipic) .and. (xip .lt. xie)) then
        fact1 = ae-apic
        fact2 = (xip-xipic)/(xie-xipic)
        axip = apic+fact1*fact2
        sxip = un-fact2
!
        fact3 = sigc/sigp1
        fact4 = (mpic/fact3+un)**(apic/axip)
!
        mxip = fact3*(fact4-sxip)
! ================================================================
! CALCUL DES VARIABLES D'ECROUISSAGES POUR LE CAS XIE< XIP < XIULT
! ================================================================
    else if ((xip .ge. xie) .and. (xip .lt. xiult)) then
        fact1 = log(un+un/eta*(xip-xie)/(xiult-xie))
        fact2 = (un-ae)/log(un+un/eta)
        axip = ae+fact1*fact2
!
        sxip = zero
!
        fact3 = sigc/sigp2
        fact4 = (me/fact3)**(ae/axip)
!
        mxip = fact3*fact4
!  ================================================================
! CALCUL DES VARIABLES D'ECROUISSAGES POUR LE CAS XIP > XIULT ----
! ================================================================
    else if (xip .ge. xiult) then
        axip = un
        sxip = zero
        mxip = mult
    end if
! ================================================================
! --- STOCKAGE ---------------------------------------------------
! ================================================================
    paraep(1) = axip
    paraep(2) = sxip
    paraep(3) = mxip
! ================================================================
end subroutine
