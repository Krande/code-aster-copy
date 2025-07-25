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
subroutine lkdndx(nmat, mater, i1, devsig, bprime, &
                  val, para, xi, dpardx, dndxi)
! person_in_charge: alexandre.foucault at edf.fr
    implicit none
!     --------------------------------------------------------------
!     CALCUL DU TERME DE LETK = DN/DXI
!     IN  NMAT     : DIMENSION TABLE DES PARAMETRES MATERIAU
!         MATER    : TABLE DES PARAMETRES MATERIAU
!         I1       : TRACE DU TENSEUR DES CONTRAINTES
!         DEVISG   : DEVIATEUR DU TENSEUR DES CONTRAINTES
!         BPRIME   : PARAMETRE DE DILATANCE FCTN SIGMA
!         VAL      : BOOLEEN SUR DILATANCE EN PRE(0) OU POST-PIC(1)
!         XI       : VARIABLE D'EXROUISSAGE XI(P OU VP)
!         PARA     : CONTIENT VALEURS DE A(XI),S(XI),M(XI)
!         DPARDX   : DERIVEE DE A(XI),S(XI),M(XI) PAR RAPPORT A XI
!     OUT DNDXI    :  DN/DXI
!     --------------------------------------------------------------
#include "asterc/r8pi.h"
#include "asterfort/cos3t.h"
#include "asterfort/lkhtet.h"
    integer(kind=8) :: nmat, val
    real(kind=8) :: i1, devsig(6), dndxi(6), bprime, mater(nmat, 2)
    real(kind=8) :: para(3), dpardx(3), xi
!
    integer(kind=8) :: i, ndt, ndi
    real(kind=8) :: sii, trois, deux, vident(6), un, zero, dbpdxi, six
    real(kind=8) :: sinpsi, dsindx, pi, pref, sigc, h0ext, s0, mult
    real(kind=8) :: xie, mvmax, mu0v, xi0v, mu1, xi1, alres, rcos3t
    real(kind=8) :: h0e, h0c, htheta, fact1, c, phi, troisd, tiers, fact2
    real(kind=8) :: sigmin, sigmax, siglim, alpha, daldxi
    real(kind=8) :: dctdxi, dtpdxi, dphidx, hx, dhxdxi, dgxdxi, dftdxi
    real(kind=8) :: daxdxi, dmxdxi, dsxdxi, dstdxi, lgleps
    real(kind=8) :: sigtil
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    parameter(six=6.d0)
    parameter(lgleps=1.0d-8)
!     --------------------------------------------------------------
    common/tdim/ndt, ndi
!     --------------------------------------------------------------
!
! --------------------------------------
! --- CONSTRUCTION VARIABLES TEMPORAIRES
! --------------------------------------
! --- VECTEUR IDENTITE
    vident(:) = zero
    do i = 1, ndi
        vident(i) = un
    end do
!
    daxdxi = dpardx(1)
    dsxdxi = dpardx(2)
    dmxdxi = dpardx(3)
!
! =================================================================
! --- CALCUL DE SII -----------------------------------------------
! =================================================================
    sii = norm2(devsig(1:ndt))
! =====================================================================
! --- RECUPERATION DE PARAMETRES DU MODELE ----------------------------
! =====================================================================
    pi = r8pi()
    pref = mater(1, 2)
    sigc = mater(3, 2)
    h0ext = mater(4, 2)
    s0 = mater(11, 2)
    mult = mater(15, 2)
    xie = mater(17, 2)
    mvmax = mater(19, 2)
!
    mu0v = mater(24, 2)
    xi0v = mater(25, 2)
    mu1 = mater(26, 2)
    xi1 = mater(27, 2)
! =================================================================
! --- CALCUL DE ALPHA RES -----------------------------------------
! =================================================================
    alres = un+mult
! =================================================================
! --- CALCUL DE H(THETA), H0E ET H0C -----------------------------
! =================================================================
    rcos3t = cos3t(devsig, pref, lgleps)
    call lkhtet(nmat, mater, rcos3t, h0e, h0c, &
                htheta)
! =================================================================
! --- CALCUL DE C TILDE -------------------------------------------
! =================================================================
    if (para(2) .le. zero) then
        fact1 = zero
        c = zero
    else
        fact1 = un+para(1)*para(3)*para(2)**(para(1)-un)
        c = sigc*(para(2))**para(1)/deux/sqrt(fact1)
    end if
! =================================================================
! --- CALCUL DE PHI TILDE -----------------------------------------
! =================================================================
    fact1 = sqrt(fact1)
    phi = deux*atan2(fact1, un)-pi/deux
! =================================================================
! --- CALCUL DE SIGMA TILDE ---------------------------------------
! =================================================================
    if (xi .le. xie) sigtil = c/tan(phi)
    if (xi .gt. xie) sigtil = zero
! =================================================================
! --- CALCUL DE SIGMIN ET SIGMAX ----------------------------------
! =================================================================
    troisd = trois/deux
    tiers = un/trois
    fact2 = (deux*htheta-(h0c+h0ext))/deux/(h0c-h0ext)
    sigmin = tiers*(i1-(troisd-fact2)*sqrt(troisd)*sii)
    sigmax = tiers*(i1+(troisd+fact2)*sqrt(troisd)*sii)
! =================================================================
! --- CALCUL DE SIGLIM  -------------------------------------------
! =================================================================
    siglim = sigmin+sigc*(mvmax*sigmin/sigc+s0)
! =================================================================
! --- CALCUL DE ALPHA  --------------------------------------------
! =================================================================
    alpha = (sigmax+sigtil)/(sigmin+sigtil)
! =================================================================
! --- CALCUL DE DSINDXI -------------------------------------------
! =================================================================
    if (val .eq. 1) then
        sinpsi = mu1*(alpha-alres)/(xi1*alpha+alres)
        if (para(2) .gt. zero) then
            dftdxi = sigc*para(2)**(para(1))*(daxdxi*log(para(2))+para(1)/para(2)*dsxdxi)
            hx = sqrt(un+para(1)*para(3)*para(2)**(para(1)-un))
            dgxdxi = daxdxi*para(3)*para(2)**(para(1)-un)+para(1)*dmxdxi*para(2)**(para(1)-un)+&
                     & para(1)*para(3)*(daxdxi*log(para(2))+(para(1)-un)/para(2)*dsxdxi)*para(2&
                     &)**(para(1)-un)
            dhxdxi = dgxdxi/(deux*hx)
            dctdxi = (dftdxi*hx-dhxdxi*sigc*para(2)**para(1))/(deux*hx**2)
            dphidx = deux/(un+hx**2)*dhxdxi
            dtpdxi = (un+(tan(phi))**2)*dphidx
            dstdxi = dctdxi/tan(phi)-c*dtpdxi/(tan(phi))**2
            daldxi = (sigmin-sigmax)/(sigmin+sigtil)**2*dstdxi
            dsindx = mu1*alres*(un+xi1)/(xi1*alpha+alres)**2*daldxi
        else
            dsindx = zero
        end if
    else
        sinpsi = mu0v*((sigmax-siglim)/(xi0v*sigmax+siglim))
        dsindx = zero
    end if
!
! --------------------------
! --- CONSTRUCTION DE DBPDXI
! --------------------------
    dbpdxi = -six*sqrt(six)/(trois-sinpsi)**2*dsindx
! -----------------------
! --- ASSEMBLAGE DE DNDXI
! -----------------------
    do i = 1, ndt
        dndxi(i) = ( &
                   devsig(i)/sii*(bprime**2+trois)-deux*bprime**2*devsig(i)/sii+deux*bprime*vide&
                   &nt(i))/(bprime**2+trois)**(trois/deux &
                   )*dbpdxi
    end do
!
end subroutine
