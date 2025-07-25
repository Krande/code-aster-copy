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
function bprime(nbmat, mater, parame, invar1, s, &
                epssig)
!
    implicit none
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/cos3t.h"
#include "asterfort/lglord.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    integer(kind=8) :: nbmat
    real(kind=8) :: mater(nbmat, 2), parame(5), invar1, s(6), epssig, bprime
! --- BUT : CALCUL DU PARAMETRE BPRIME ---------------------------------
! ======================================================================
! IN  : NBMAT  : NOMBRE DE PARAMETRES DU MODELE ------------------------
! --- : MATER  : PARAMETRES DU MODELE ----------------------------------
! --- : PARAME : VARIABLES D'ECROUISSAGE -------------------------------
! --- : NDT    : NOMBRE DE COMPOSANTES TOTALES DU TENSEUR --------------
! --- : INVAR1 : PREMIER INVARIANT DU TENSEUR DES CONTRAINTES ----------
! --- : S      : DEVIATEUR DES CONTRAINTES -----------------------------
! --- : EPSSIG : EPSILON -----------------------------------------------
! OUT : BPRIME : PARAMETRE CONTROLANT LE COMPORTEMENT VOLUMIQUE --------
! ------------ : DU MATERIAU -------------------------------------------
! ======================================================================
    integer(kind=8) :: ndt, ndi
    real(kind=8) :: mun, un, deux, trois, six, epstol, un_m_rcos2
    real(kind=8) :: mult, sigc, gamma, ksi, pref, prec
    real(kind=8) :: sgamp, agamp, mgamp, sii, fact1, fact2
    real(kind=8) :: rcos3t
    real(kind=8) :: phi0, c0, sigt0, sig1, sig2, sig3, alpha, sinpsi
    blas_int :: b_incx, b_incy, b_n
! ======================================================================
! --- INITIALISATION DE PARAMETRES -------------------------------------
! ======================================================================
    parameter(mun=-1.0d0)
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
    parameter(six=6.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
    epstol = r8prem()
    sigt0 = 0.0d0
! ======================================================================
! --- INITIALISATION DES PARAMETRES MATERIAU ---------------------------
! ======================================================================
    mult = mater(3, 2)
    sigc = mater(9, 2)
    gamma = mater(10, 2)
    ksi = mater(11, 2)
    pref = mater(15, 2)
! ======================================================================
! --- RECUPERATION DES VARIABLES D'ECROUISSAGE -------------------------
! ======================================================================
    sgamp = parame(1)
    agamp = parame(2)
    mgamp = parame(4)
! ======================================================================
! --- CALCULS INTERMEDIAIRE POUR LE CALCUL DE BPRIME -------------------
! ======================================================================
    b_n = to_blas_int(ndt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    sii = ddot(b_n, s, b_incx, s, b_incy)
    sii = sqrt(sii)
    rcos3t = cos3t(s, pref, epssig)
! ======================================================================
! --- CALCUL DE PHI0 = 2*ARCTAN(RAC(1+A*M*S**(A-1))) - PI/2 ------------
! ======================================================================
! --- SI SGAMP = 0 ON A ALPHA = SIG1/SIG3 ------------------------------
! ======================================================================
    if (sgamp .lt. epstol) goto 10
    fact1 = sgamp**(agamp-un)
    fact2 = un+agamp*mgamp*fact1
    if (fact2 .lt. epstol) then
        call utmess('F', 'ALGELINE_4')
    end if
    fact2 = sqrt(fact2)
    phi0 = deux*atan2(fact2, un)-r8pi()/deux
! ======================================================================
! --- CALCUL DE C0 = SIGC*S**A/(RAC(1+A*M*S**(A-1))) -------------------
! ======================================================================
    c0 = sigc*sgamp**agamp/fact2
! ======================================================================
! --- CALCUL DE SIGT0 = 2*C0*RAC((1-SIN(PHI0))/(1+SIN(PHI0)) -----------
! ======================================================================
    if ((un+sin(phi0)) .lt. epstol) then
        call utmess('F', 'ALGELINE_4')
    end if
    sigt0 = deux*c0*sqrt((un-sin(phi0))/(un+sin(phi0)))
10  continue
! ======================================================================
! --- CALCULS DE INTERMEDIAIRE -----------------------------------------
! ======================================================================
    sig1 = invar1/trois+sqrt(deux/trois)*sii*rcos3t
!
!   Variabilite machine: on met a 0 si un-rcos3t*rcos3t trop petit
    un_m_rcos2 = un-rcos3t*rcos3t
    prec = 100.d0*r8prem()
!
    if (un_m_rcos2 .lt. prec) then
        un_m_rcos2 = 0.d0
    end if
!
    sig2 = invar1/trois-sqrt(deux/trois)*sii*(rcos3t/deux+sqrt(trois*(un_m_rcos2))/deux)
    sig3 = invar1/trois+sqrt(deux/trois)*sii*(-rcos3t/deux+sqrt(trois*(un_m_rcos2))/deux)
! ======================================================================
! --- RECUPERATION DE SIG1 (MAX) ET SIG3 (MIN) -------------------------
! ======================================================================
    call lglord(sig1, sig2, sig3)
! ======================================================================
! --- CALCUL DE ALPHA = (SIG1-SIGT0)/(SIG3-SIGT0) ----------------------
! ======================================================================
    if (abs(sig3-sigt0) .lt. epstol) then
        alpha = un
    else
        alpha = (sig1-sigt0)/(sig3-sigt0)
    end if
! ======================================================================
! --- CALCUL DE SIN(PSI) = GAMMA*(ALPHA-MULT-1) / (KSI*ALPHA+MULT+1) ---
! ======================================================================
    sinpsi = gamma*(alpha-mult-un)/(ksi*alpha+mult+un)
! ======================================================================
! --- AJUSTEMENT DE LA LOI DE DILATANCE --------------------------------
! ======================================================================
    if (sinpsi .lt. mun) sinpsi = mun
    if (sinpsi .gt. un) sinpsi = un
! ======================================================================
! --- CALCUL DE BPRIME = -2*RAC(6)*SIN(PSI)/(3-SIN(PSI)) ---------------
! ======================================================================
    bprime = mun*deux*sqrt(six)*sinpsi/(trois-sinpsi)
! ======================================================================
end function
