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

subroutine lcmmfe(taus, coeft, materf, ifa, nmat, &
                  nbcomm, necoul, is, nbsys, vind, &
                  dy, rp, alphap, gammap, dt, &
                  dalpha, dgamma, dp, crit, sgns, &
                  nfs, nsg, hsr, iret)
! aslint: disable=W1504
    implicit none
#include "asterc/r8miem.h"
#include "asterfort/lcddcc.h"
#include "asterfort/lcmmdd.h"
#include "asterfort/lcmmkr.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ifa, nmat, nbcomm(nmat, 3), iret
    integer(kind=8) :: ifl, is, nbsys, nuecou, nfs, nsg
    real(kind=8) :: taus, coeft(nmat), alphap, dgamma, dp, dt
    real(kind=8) :: rp, sgns, hsr(nsg, nsg), dy(*), vind(*), materf(nmat)
    real(kind=8) :: dalpha
    real(kind=8) :: c, k, n, ftau, crit, a, d
    real(kind=8) :: gammap, ptit, cisa2
    character(len=16) :: necoul
    integer(kind=8) :: irr, decirr, nums, decal, gdef
    common/polycr/irr, decirr, nums, decal, gdef
! person_in_charge: jean-michel.proix at edf.fr
! ======================================================================
!  COMPORTEMENT MONOCRISTALLIN : ECOULEMENT (VISCO)PLASTIQUE
!  INTEGRATION DES LOIS MONOCRISTALLINES
!       IN  TAUS    :  SCISSION REDUITE
!           COEFT   :  PARAMETRES MATERIAU
!           IFA     :  NUMERO DE FAMILLE
!           NMAT    :  NOMBRE MAXI DE MATERIAUX
!           NBCOMM  :  NOMBRE DE COEF MATERIAU PAR FAMILLE
!           NECOUL  :  NOM DE LA LOI D'ECOULEMENT
!           RP      :  R(P) FONCTION D'ECROUISSAGE ISTROPE
!           ALPHAP  :  ALPHA A T ACTUEL
!           GAMMAP  :  GAMMA A T ACTUEL
!           DT      :  INTERVALLE DE TEMPS EVENTULLEMENT REDECOUPE
!     OUT:
!           DGAMMA  :  DEF PLAS
!           DALPHA  :  VARIABLE dalpha pour Kocks-Rauch
!           DP      :  DEF PLAS CUMULEE
!           CRIT    :  CRITERE
!           SGNS    :  SIGNE DE GAMMA
!           IRET    :  CODE RETOUR
! ======================================================================
!
!
!     DANS VIS : 1 = ALPHA, 2=GAMMA, 3=P
!
    ifl = nbcomm(ifa, 1)
    nuecou = nint(coeft(ifl))
    ptit = r8miem()
    dgamma = 0.d0
    dp = 0.d0
    crit = 0.d0
    sgns = 0.d0
    iret = 0
!
!-------------------------------------------------------------
!     POUR UN NOUVEAU TYPE D'ECOULEMENT, CREER UN BLOC IF
!------------------------------------------------------------
!
!      IF (NECOUL.EQ.'MONO_VISC1') THEN
    if (nuecou .eq. 1) then
        n = coeft(ifl+1)
        k = coeft(ifl+2)
        c = coeft(ifl+3)
        ftau = taus-c*alphap
        if (abs(ftau) .lt. ptit) then
            sgns = 1.d0
        else
            sgns = ftau/abs(ftau)
        end if
        crit = abs(ftau)-rp
        if (crit .gt. 0.d0) then
            dp = ((crit/k)**n)*dt
            dgamma = dp*sgns
        else
            dp = 0.d0
            dgamma = 0.d0
        end if
!
!      IF (NECOUL.EQ.'MONO_VISC2') THEN
    else if (nuecou .eq. 2) then
        n = coeft(ifl+1)
        k = coeft(ifl+2)
        c = coeft(ifl+3)
        a = coeft(ifl+4)
        d = coeft(ifl+5)
!
        ftau = taus-c*alphap-a*gammap
        if (abs(ftau) .lt. ptit) then
            sgns = 1.d0
        else
            sgns = ftau/abs(ftau)
        end if
!
        crit = abs(ftau)-rp+0.5d0*d*c*alphap**2
        if (crit .gt. 0.d0) then
            dp = ((crit/k)**n)*dt
            dgamma = dp*sgns
        else
            dp = 0.d0
            dgamma = 0.d0
        end if
!
!      IF (NECOUL.EQ.'KOCKS_RAUCH') THEN
    else if (nuecou .eq. 4) then
        if (materf(nmat) .eq. 0) then
            cisa2 = (materf(1)/2.d0/(1.d0+materf(2)))**2
        else
            cisa2 = (materf(36)/2.d0)**2
        end if
        call lcmmkr(taus, coeft, cisa2, ifa, nmat, &
                    nbcomm, is, nbsys, nfs, nsg, &
                    hsr, vind(decal+1), dy, dt, dalpha, &
                    dgamma, dp, crit, sgns, iret)
!         CALCUL D'UN RP FICTIF POUR LCMMVX :
!         CELA PERMET D'ESTIMER LE PREMIER POINT DE NON LINEARITE
        rp = coeft(ifl+3)
!
!      IF (NECOUL.EQ.'MONO_DD_CFC') THEN + IRRA
    else if ((nuecou .eq. 5) .or. (nuecou .eq. 6) .or. (nuecou .eq. 8)) then
        call lcmmdd(taus, coeft, ifa, nmat, nbcomm, &
                    is, nbsys, nfs, nsg, hsr, &
                    vind(decal+1), dy, dt, rp, nuecou, &
                    dalpha, dgamma, dp, iret)
        if (iret .gt. 0) goto 9999
        crit = dp
        if (abs(taus) .le. ptit) then
            sgns = 1.d0
        else
            sgns = taus/abs(taus)
        end if
!
!      IF (NECOUL.EQ.'MONO_DD_CC')OR (NECOUL.EQ.'MONO_DD_CC_IRRA') THEN
    else if (nuecou .eq. 7) then
        call lcddcc(taus, coeft, ifa, nmat, nbcomm, &
                    is, nbsys, nfs, nsg, hsr, &
                    vind, dy, dt, rp, nuecou, &
                    dalpha, dgamma, dp, iret)
        if (iret .gt. 0) goto 9999
        crit = dp
        if (abs(taus) .le. ptit) then
            sgns = 1.d0
        else
            sgns = taus/abs(taus)
        end if
!
!
!
!
    else
        call utmess('F', 'COMPOR1_20')
    end if
9999 continue
end subroutine
