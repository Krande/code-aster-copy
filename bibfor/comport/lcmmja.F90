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

subroutine lcmmja(typmod, nmat, materf, timed, &
                  timef, itmax, toler, nbcomm, cpmono, &
                  pgl, nfs, nsg, toutms, hsr, &
                  nr, nvi, vind, df, yf, &
                  yd, dy, drdy, iret)
! aslint: disable=W1306,W1504
    implicit none
!       MONOCRISTAL : CALCUL DU JACOBIEN DU SYSTEME NL A RESOUDRE = DRDY
!                    DY    = ( DSIG + DGAMMA PAR SYST )
!                    Y     = ( SIG   GAMMA P par syst. gliss)
!       IN  COMP   :  NOM COMPORTEMENT
!           TYPMOD :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           TIMED  :  ISTANT PRECEDENT
!           TIMEF  :  INSTANT ACTUEL
!           ITMAX  :  ITER_INTE_MAXI
!           TOLER  :  RESI_INTE
!           NBCOMM :  INCIDES DES COEF MATERIAU
!           CPMONO :  NOM DES COMPORTEMENTS
!           PGL    :  MATRICE DE PASSAGE
!           TOUTMS :  TENSEURS D'ORIENTATION
!           HSR    :  MATRICE D'INTERACTION
!           NR     :  DIMENSION DECLAREE DRDY
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!           DF     :  Increment de Gradient de deformation
!           YD     :  VARIABLES A T
!           YF     :  VARIABLES A T + DT
!           DY     :  SOLUTION
!       OUT DRDY   :  JACOBIEN DU SYSTEME NON LINEAIRE
!           IRET   :  CODE RETOUR
!       ----------------------------------------------------------------
#include "asterfort/calcfe.h"
#include "asterfort/caldfe.h"
#include "asterfort/caldto.h"
#include "asterfort/caltau.h"
#include "asterfort/lcicma.h"
#include "asterfort/lcmmjb.h"
#include "asterfort/lcmmjg.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/lcopil.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: nmat, nr, nbfsys, ndt, ndi, nsfa, nsfv, nbsys, is, ir
    integer(kind=8) :: nbcomm(nmat, 3), ifa, i, j, k, l, iret, ifl, itmax, nuvr, nuvs
    integer(kind=8) :: nuecou, ind(3, 3), nvi, nfs, nsg, iexp
    real(kind=8) :: vind(*), yf(*), dy(*), drdy(nr, nr), materf(nmat*2)
    real(kind=8) :: pgl(3, 3), toutms(nfs, nsg, 6), hsr(nsg, nsg), gamsns(3, 3)
    real(kind=8) :: timed, timef, msdgdt(6, 6), dt, fkooh(6, 6), sigf(6)
    real(kind=8) :: toler, dgsdts, dksdts, dgrdbs, dkrdbs, taus, taur
    real(kind=8) :: msns(3, 3)
    real(kind=8) :: q(3, 3), mus(6), ns(3), ms(3), mur(6), dtods(3, 3)
    real(kind=8) :: dfpds(3, 3, 3, 3), yd(*), msnst(3, 3, nsg), fp(3, 3)
    real(kind=8) :: mrnr(3, 3), df(3, 3), fe(3, 3), expbp(nsg)
    real(kind=8) :: dfpdbs(3, 3, nsg), dfpdga(3, 3, nsg)
    character(len=16) :: nomfam
    character(len=24) :: cpmono(5*nmat+1)
    character(len=8) :: typmod
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
    data ind/1, 4, 5, 4, 2, 6, 5, 6, 3/
!     ----------------------------------------------------------------
!
    iret = 0
    dt = timef-timed
!
    call r8inir(nr*nr, 0.d0, drdy, 1)
    call r8inir(36, 0.d0, msdgdt, 1)
!
    sigf(1:ndt) = yf(1:ndt)
!
!     Inverse de la matrice de Hooke
    if (materf(nmat) .eq. 0) then
        call lcopil('ISOTROPE', typmod, materf(1), fkooh)
    else if (materf(nmat) .eq. 1) then
        call lcopil('ORTHOTRO', typmod, materf(1), fkooh)
    end if
!
    if (gdef .eq. 1) then
        call r8inir(81, 0.d0, dfpds, 1)
        call r8inir(3*3*nsg, 0.d0, dfpdbs, 1)
!        calcul de DFPDGA : dFp / dGamma_S pour tous les systemes S
        call lcmmjg(nmat, nbcomm, cpmono, hsr, &
                    dt, nvi, vind, yd, dy, &
                    itmax, toler, materf, sigf, fkooh, &
                    nfs, nsg, toutms, pgl, msnst, &
                    gamsns, dfpdga, iret)
    end if
!
!     NSFA : debut de la famille IFA dans DY et YD, YF
    nsfa = 6
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
!     LE NUMERO GLOBAL DU SYSTEME IS DANS Y EST NUVS
    nbfsys = nbcomm(nmat, 2)
!
    do ifa = 1, nbfsys
!
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        ifl = nbcomm(ifa, 1)
        nuecou = nint(materf(nmat+ifl))
!
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ns, ms, 0, q)
!
        do is = 1, nbsys
!
!           calcul de Tau_s HPP ou GDEF
!
            call caltau(ifa, is, sigf, fkooh, &
                        nfs, nsg, toutms, taus, mus, &
                        msns)
!
            nuvs = nsfa+is
!
!           CALCUL DES DERIVEES :
!           DGSDTS=dGamma_S/dTau_S,  DKSDTS=dK_s/dTau_S,
!           DGRDBS=dGamma_R/dBeta_S, DKRDBS=dK_S/dBeta_R
!
            iexp = 0
            if (is .eq. 1) iexp = 1
            call lcmmjb(taus, materf, cpmono, ifa, nmat, &
                        nbcomm, dt, nuecou, nsfv, nsfa, &
                        is, is, nbsys, nfs, nsg, &
                        hsr, vind, dy, iexp, expbp, &
                        itmax, toler, dgsdts, dksdts, dgrdbs, &
                        dkrdbs, iret)
!           ici  DGRDBS,DKRDBS sont inutiles
            if (iret .gt. 0) goto 999
!
            if (abs(dgsdts) .gt. 0.d0) then
!              Cas ou Delta-Gamma_S est non nul
                if (gdef .eq. 0) then
!                 dR1/dS
                    do i = 1, 6
                        do j = 1, 6
                            msdgdt(i, j) = msdgdt(i, j)+mus(i)*mus(j)* &
                                           dgsdts
                        end do
                    end do
!                 dR2/dS
                    do i = 1, 6
                        drdy(nuvs, i) = -mus(i)*dksdts
                    end do
                else
                    call caldto(sigf, fkooh, msns, dtods)
!                 dR1/dS
                    do i = 1, 3
                        do j = 1, 3
                            do k = 1, 3
                                do l = 1, 3
                                    dfpds(i, j, k, l) = dfpds(i, j, k, l)+ &
                                                        dfpdga(i, j, is)*dgsdts*dtods(k, l)
                                end do
                            end do
                        end do
                    end do
!                 dR2/dS
                    do i = 1, 3
                        do j = 1, 3
                            drdy(nuvs, ind(i, j)) = -dksdts*dtods(i, j)
                        end do
                    end do
                end if
            end if
!
!------------------------
!           calcul des ns termes dR1_i/dBeta_s
!           et     des ns termes dR2_r/dBeta_s
!------------------------
            do ir = 1, nbsys
                call caltau(ifa, ir, sigf, fkooh, &
                            nfs, nsg, toutms, taur, mur, &
                            mrnr)
!
                nuvr = nsfa+ir
!
                call lcmmjb(taur, materf, cpmono, ifa, nmat, &
                            nbcomm, dt, nuecou, nsfv, nsfa, &
                            ir, is, nbsys, nfs, nsg, &
                            hsr, vind, dy, iexp, expbp, &
                            itmax, toler, dgsdts, dksdts, dgrdbs, &
                            dkrdbs, iret)
!              ici DGSDTS,DKSDTS sont inutiles
                if (iret .gt. 0) goto 999
!
                if (abs(dgrdbs) .gt. 0.d0) then
                    if (gdef .eq. 0) then
!                    terme dR1/dAlpha_s
                        do i = 1, 6
                            drdy(i, nuvs) = drdy(i, nuvs)+mur(i)*dgrdbs
                        end do
                    else
                        do i = 1, 3
                            do j = 1, 3
                                dfpdbs(i, j, is) = dfpdbs(i, j, is)+ &
                                                   dfpdga(i, j, ir)*dgrdbs
                            end do
                        end do
                    end if
!                 terme dR2r/dGammas
                    drdy(nuvr, nuvs) = -dkrdbs
                end if
            end do
            drdy(nuvs, nuvs) = drdy(nuvs, nuvs)+1.d0
        end do
        nsfa = nsfa+nbsys
        nsfv = nsfv+nbsys*3
    end do
!
    if (gdef .eq. 1) then
        call calcfe(nr, ndt, nvi, vind, df, &
                    gamsns, fe, fp, iret)
        call caldfe(df, nr, nvi, vind, dfpds, &
                    fe, dfpdbs, msdgdt, drdy)
    end if
!
    msdgdt(1:ndt, 1:ndt) = msdgdt(1:ndt, 1:ndt)+fkooh(1:ndt, 1:ndt)
    call lcicma(msdgdt, 6, 6, ndt, ndt, &
                1, 1, drdy, nr, nr, &
                1, 1)
999 continue
end subroutine
