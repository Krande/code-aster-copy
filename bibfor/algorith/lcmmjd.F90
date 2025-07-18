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
subroutine lcmmjd(taur, materf, ifa, nmat, nbcomm, &
                  dt, ir, is, nbsys, nfs, &
                  nsg, hsr, vind, dy, dpdtau, &
                  dprdas, dhrdas, hr, dpr, sgnr, &
                  iret)
! aslint: disable=W1306,W1504
    implicit none
#include "asterfort/assert.h"
#include "asterfort/lcmmdc.h"
#include "asterfort/lcmmdh.h"
#include "asterfort/lcmmfe.h"
#include "asterfort/lcmmfi.h"
    integer(kind=8) :: ifa, nmat, nbcomm(nmat, 3), nfs, nsg
    real(kind=8) :: taur, materf(nmat*2), rr, dt, vind(36), dy(12)
    real(kind=8) :: dpdtau, dprdas, hsr(nsg, nsg), hr
! person_in_charge: jean-michel.proix at edf.fr
! ======================================================================
!  CALCUL DES DERIVEES DES VARIABLES INTERNES DES LOIS MONOCRISTALLINES
!  POUR LA LOI D'ECOULEMENT  DD-CFC
!       IN  TAUR    :  SCISSION REDUITE
!           MATERF  :  PARAMETRES MATERIAU
!           IFA     :  NUMERO DE FAMILLE
!           NBCOMM  :  NOMBRE DE COEF MATERIAU PAR FAMILLE
!           DT      :  INCREMENT DE TEMPS
!           IR,IS   :  NUMEROS DE SYSTEMES DE GLISSEMENT
!           NBSYS   :  NOMBRE DE SYSTEMES DE GLISSEMENT (12 POUR CFC)
!           HSR     :  MATRICE D'INTERACTION
!           VIND,DY :  VARIABLES INTERNES A T ET SOLUTION A T+DT
!     OUT:
!           DPDTAU  :  dpr/dTaur
!           DPRDAS  :  dpr/dAlphas
!           DHRDAS  :  dHr/dAlphaS
!     ----------------------------------------------------------------
    real(kind=8) :: b, n, a, gamma0, r8b, dpr, dtrdas, critr, expbp(nsg)
    real(kind=8) :: ceff, dcdals, soms3, soms2, soms1
    real(kind=8) :: dhrdas, tauf, y, beta, mu, somaal, ars, sgnr, t3
    real(kind=8) :: alphap(12)
    integer(kind=8) :: ifl, is, nbsys, ir, iret, nuecou, iei, iu, i, is3, ir3
    character(len=16) :: k16b
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
!
    ifl = nbcomm(ifa, 1)+nmat
    nuecou = nint(materf(ifl))
    iei = nbcomm(ifa, 3)+nmat
    iret = 0
!
    if ((nuecou .ne. 5) .and. (nuecou .ne. 8)) then
        ASSERT(.false.)
    end if
!             MATRICE JACOBIENNE DU SYSTEME :
!  R1 = D-1*SIGMA - (D_M-1*SIGMA_M)-(DEPS-DEPS_TH)+Somme(ms*Dps*S)=0
!  R2 = dALPHA - Dps*h(alphas)
! avec S=sgn(TAUR)
!
! ON VEUT CALCULER :
!     d(Dps)/dTAUR
    dpdtau = 0.d0
!     d(R1)/dalphas
    hr = 0.d0
!     d(R2)/d(alphar)
    dhrdas = 0.d0
!     DPSDAR=d(Dp_s)/d(Alpha_r)
    dprdas = 0.d0
!
    tauf = materf(ifl+1)
    gamma0 = materf(ifl+2)
    a = materf(ifl+3)
    b = materf(ifl+4)
    n = materf(ifl+5)
    y = materf(ifl+6)
    beta = materf(iei+2)
    if (nuecou .eq. 5) then
        mu = materf(iei+4)
    else if (nuecou .eq. 8) then
        mu = materf(iei+12)
    end if
    k16b = ' '
!     CALCUL de l'écrouissage RR=TAUr_Forest
    call lcmmfi(materf(nmat+1), ifa, nmat, nbcomm, k16b, &
                ir, nbsys, vind, decal, dy, &
                nfs, nsg, hsr, 1, expbp, &
                rr)
    if (iret .gt. 0) goto 999
!
!     CALCUL de l'écoulement dpr et du critère
    call lcmmfe(taur, materf(nmat+1), materf(1), ifa, nmat, &
                nbcomm, k16b, ir, nbsys, vind, &
                dy, rr, r8b, r8b, dt, &
                r8b, r8b, dpr, critr, sgnr, &
                nfs, nsg, hsr, iret)
    if (iret .gt. 0) goto 999
!
!
!     1. d(Dp_r)/d(Tau_s)
    if (dpr .gt. 0.d0) then
        if (abs(taur) .gt. 0.d0) then
            dpdtau = n*(dpr+gamma0*dt)/taur
        end if
    end if
!
    do iu = 1, nbsys
        alphap(iu) = vind(decal+3*(iu-1)+1)+dy(iu)
    end do
!
    call lcmmdc(materf(nmat+1), ifa, nmat, nbcomm, alphap, &
                is, ceff, dcdals)
!
    call lcmmdh(materf(nmat+1), ifa, nmat, nbcomm, alphap, &
                nfs, nsg, hsr, nbsys, ir, &
                nuecou, hr, soms1, soms2, soms3)
!
!
    somaal = 0.d0
    do i = 1, 12
        if (alphap(i) .gt. 0.d0) then
            somaal = somaal+hsr(ir, i)*alphap(i)
        end if
    end do
!   rr ne doit pas etre nul car ce sont des densites de dislocations
    if (abs(rr) .gt. 1.e-20) then
        dtrdas = mu*mu*ceff/2.0d0/rr*(2.0*dcdals*somaal+ceff*hsr(ir, is))
    else
        iret = 1
        goto 999
    end if
!
!
!     2. d(Dp_r)/d(Omega_s)
!
    if (dpr .gt. 0.d0) then
        dprdas = -n*(dpr+gamma0*dt)/(tauf+rr)*dtrdas
    end if
!
    ars = hsr(ir, is)
    if (ars*alphap(is) .gt. 0.d0) then
        t3 = ars/2.d0/sqrt(ars*alphap(is))
    end if
!
    dhrdas = 0.d0
!     IS APPARTIENT-IL A FOREST(IR) ?
!     division entiere
    is3 = (is-1)/3
    ir3 = (ir-1)/3
    if (is3 .ne. ir3) then
        if (alphap(is) .gt. 0.d0) then
            dhrdas = dhrdas+a*(sqrt(ars)/soms1)
        end if
    end if
!
    if (soms1 .gt. 0.d0) then
        dhrdas = dhrdas-a*t3*soms2/soms1/soms1
    end if
!     IS APPARTIENT-IL A COPLA(IR) ?
    if (is3 .eq. ir3) then
        if (ars*alphap(is) .gt. 0.d0) then
            dhrdas = dhrdas+b*t3*ceff
        end if
    end if
!
    dhrdas = dhrdas+b*soms3*dcdals
!
!     3. d(h_r)/d(Omega_s)
    if (is .eq. ir) dhrdas = dhrdas-y/beta
!
999 continue
!
end subroutine
