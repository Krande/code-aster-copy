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

subroutine lcmmj1(taur, materf, cpmono, ifa, nmat, &
                  nbcomm, dt, nsfv, nsfa, ir, &
                  is, nbsys, nfs, nsg, hsr, &
                  vind, dy, iexp, expbp, itmax, &
                  toler, dgsdts, dksdts, dgrdbs, dkrdbs, &
                  iret)
! aslint: disable=,W1504
    implicit none
! person_in_charge: jean-michel.proix at edf.fr
!       ----------------------------------------------------------------
!       MONOCRISTAL : DERIVEES DES TERMES UTILES POUR LE CALCUL
!                    DU JACOBIEN DU SYSTEME NL A RESOUDRE = DRDY
!                    cf. R5.03.11 comportements MONO_VISC*
!       IN
!           TAUR   :  SCISSION REDUITE SYSTEME IR
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           CPMONO :  NOM DES COMPORTEMENTS
!           IFA    :  NUMERO FAMILLE
!           NMAT   :  DIMENSION MATER
!           NBCOMM :  INCIDES DES COEF MATERIAU
!           DT     :  ACCROISSEMENT INSTANT ACTUEL
!           NSFV   :  DEBUT DES SYST. GLIS. DE LA FAMILLE IFA DANS VIND
!           NSFA   :  DEBUT DES SYST. GLIS. DE LA FAMILLE IFA DANS Y
!           IS     :  NUMERO DU SYST. GLIS. S
!           IR     :  NUMERO DU SYST. GLIS. R
!           NBSYS  :  NOMBRE DE SYSTEMES DE GLISSEMENT FAMILLE IFA
!           HSR    :  MATRICE D'INTERACTION
!           VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!           DY     :  SOLUTION           =  ( DSIG DX1 DX2 DP (DEPS3) )
!           ITMAX  :  ITER_INTE_MAXI
!           TOLER  :  RESI_INTE
!       OUT DGSDTS :  derivee dGammaS/dTauS
!       OUT DKSDTS :  dkS/dTaus
!       OUT DGRDBS :  dGammaR/dBetaS
!       OUT DKRDBS :  dkR/dBetaS
!       OUT IRET   :  CODE RETOUR
!       ----------------------------------------------------------------
#include "asterfort/lcmmfc.h"
#include "asterfort/lcmmfe.h"
#include "asterfort/lcmmfi.h"
#include "asterfort/lcmmjc.h"
#include "asterfort/lcmmjf.h"
#include "asterfort/lcmmji.h"
    integer(kind=8) :: nmat, nuvr, nuvs, iexp, ir, nsfa, nsfv, itmax, nfs, nsg
    integer(kind=8) :: nbcomm(nmat, 3), nuvi, ifa, nbsys, is, iret
    real(kind=8) :: vind(*), dgdtau, dgrdrr
    real(kind=8) :: hsr(nsg, nsg), expbp(nsg), dp, dy(*)
    real(kind=8) :: materf(nmat*2), dt, dgamms, rr
    real(kind=8) :: alpham, dalpha, alphar, crit, dardgr, dgamma
    real(kind=8) :: sgns, taur, gammar, pms
    real(kind=8) :: dgrdar, toler, drrdps, ps, petith
    real(kind=8) :: dgsdts, dksdts, dgrdbs, dkrdbs
    character(len=24) :: cpmono(5*nmat+1)
    character(len=16) :: necoul, necris, necrci
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
    iret = 0
    dgsdts = 0.d0
    dksdts = 0.d0
    dgrdbs = 0.d0
    dkrdbs = 0.d0
!
    necoul = cpmono(5*(ifa-1)+3) (1:16)
    necris = cpmono(5*(ifa-1)+4) (1:16)
    necrci = cpmono(5*(ifa-1)+5) (1:16)
!
    nuvr = nsfa+ir
    nuvs = nsfa+is
    nuvi = nsfv+3*(ir-1)
!      PM=VIND(NUVI+3)
    alpham = vind(nuvi+1)
    dgamms = dy(nuvr)
    gammar = vind(nuvi+2)+dgamms
!
!     CALCUL DE DALPHA
    call lcmmfc(materf(nmat+1), ifa, nmat, nbcomm, necrci, &
                itmax, toler, alpham, dgamms, dalpha, &
                iret)
    if (iret .gt. 0) goto 9999
    alphar = alpham+dalpha
!
!     CALCUL DE R(P) : RP=R0+Q*(1.D0-EXP(-B*P))
!        ECROUISSAGE ISOTROPE : CALCUL DE R(P)
!          IEXP=1
!          IF (IR.EQ.1) IEXP=1
    call lcmmfi(materf(nmat+1), ifa, nmat, nbcomm, necris, &
                ir, nbsys, vind, nsfv, dy(nsfa+1), &
                nfs, nsg, hsr, iexp, expbp, &
                rr)
!
!     CALCUL de DGAMMA et de CRIT
    decal = nsfv
!
    call lcmmfe(taur, materf(nmat+1), materf(1), ifa, nmat, &
                nbcomm, necoul, ir, nbsys, vind, &
                dy(nsfa+1), rr, alphar, gammar, dt, &
                dalpha, dgamma, dp, crit, sgns, &
                nfs, nsg, hsr, iret)
    if (iret .gt. 0) goto 9999
!
    if (crit .gt. 0.d0) then
!        CALCUL de dF/dtau
        call lcmmjf(taur, materf(nmat+1), materf(1), ifa, nmat, &
                    nbcomm, dt, necoul, ir, is, &
                    nbsys, vind(nsfv+1), dy(nsfa+1), nfs, nsg, &
                    hsr, rr, alphar, dalpha, gammar, &
                    dgamms, sgns, dgdtau, dgrdar, dgrdrr, &
                    petith, iret)
        if (iret .gt. 0) goto 9999
!
        dgsdts = dgdtau
        dksdts = dgdtau
!
!        CALCUL DE dRr/dps
!         PS=PM+ABS(DY(NUVS))
        pms = vind(nsfv+3*(is-1)+3)
        ps = pms+abs(dy(nuvs))
        call lcmmji(materf(nmat+1), ifa, nmat, nbcomm, necris, &
                    nfs, nsg, hsr, ir, is, &
                    ps, drrdps)
!
!        CALCUL DE DALPHAs/dGAMMAs
        call lcmmjc(materf(nmat+1), ifa, nmat, nbcomm, ir, &
                    is, necrci, dgamms, alpham, dalpha, &
                    sgns, dardgr)
!
        dgrdbs = dgrdar*dardgr+dgrdrr*drrdps*sgns
!
        dkrdbs = dgrdbs
!
    end if
!
9999 continue
!
end subroutine
