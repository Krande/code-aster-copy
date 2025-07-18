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
subroutine lcmmjg(nmat, nbcomm, cpmono, hsr, dt, &
                  nvi, vind, yd, dy, itmax, &
                  toler, materf, sigf, fkooh, nfs, &
                  nsg, toutms, pgl, msnst, gamsns, &
                  dfpdga, iret)
! aslint: disable=W1306,W1504
    implicit none
!     MONOCRISTAL : POUR LE CALCUL DU JACOBIEN DU SYSTEME NL A RESOUDRE
!                   CALCUL EN GDEF DE Gamma.Ms*Ns et dFp/Dgamma
!       IN
!           COMP   :  NOM COMPORTEMENT
!           NMAT   :  DIMENSION MATER
!           NBCOMM :  INCIDES DES COEF MATERIAU
!           CPMONO :  NOM DES COMPORTEMENTS
!           HSR    :  MATRICE D'INTERACTION
!           DT     :  DELTA T
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!           YD     :  VARIABLES A T
!           DY     :  SOLUTION
!           ITMAX  :  ITER_INTE_MAXI
!           TOLER  :  RESI_INTE
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           SIGF   :  CONTRAINTES A T+DT
!           FKOOH  :  INVERSE TENSEUR HOOKE
!           PGL    :  MATRICE DE PASSAGE
!           TOUTMS :  TENSEURS D'ORIENTATION
!       OUT MSNST  :  Ms*Ns pour chaque systele de glissement
!           GAMSNS :  Somme de GammaS*MS*NS
!           DFPDGA :  derivee de Fp / dGamma_S pour tous les systemes S
!       OUT IRET   :  CODE RETOUR
!       ----------------------------------------------------------------
#include "asterfort/assert.h"
#include "asterfort/caldfp.h"
#include "asterfort/caltau.h"
#include "asterfort/lcmmlc.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
    integer(kind=8) :: nvi, nmat, nbfsys, nsfa, nsfv, nbsys, is, nfs, nsg
    integer(kind=8) :: nbcomm(nmat, 3), ifa, iret, itmax
    real(kind=8) :: vind(*), dy(*), materf(nmat*2)
    real(kind=8) :: pgl(3, 3), toutms(nfs, nsg, 6), hsr(nsg, nsg), gamsns(3, 3)
    real(kind=8) :: dt, fkooh(6, 6), sigf(6), toler, taus, dp, crit, sgns, rp
    real(kind=8) :: q(3, 3), mus(6), ns(3), ms(3), dfpdga(3, 3, nsg)
    real(kind=8) :: expbp(nsg), yd(*), msnst(3, 3, nsg), dalpha, dgamma
    character(len=16) :: nomfam
    character(len=24) :: cpmono(5*nmat+1)
    blas_int :: b_incx, b_incy, b_n
!     ----------------------------------------------------------------
!
!     NSFA : debut de la famille IFA dans DY et YD, YF
    nsfa = 6
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
    nbfsys = nbcomm(nmat, 2)
!     PROGRAMMATION VALABLE POUR UNE SEULE FAMILLE DE SYSTEMES
    ASSERT(nbfsys .eq. 1)
    do ifa = 1, nbfsys
!        Calcul preliminaire de somme(dgamma*ms*ns)
        call r8inir(9, 0.d0, gamsns, 1)
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ns, ms, 0, q)
        do is = 1, nbsys
            call caltau(ifa, is, sigf, fkooh, nfs, &
                        nsg, toutms, taus, mus, msnst(1, 1, is))
            call lcmmlc(nmat, nbcomm, cpmono, nfs, nsg, &
                        hsr, nsfv, nsfa, ifa, nbsys, &
                        is, dt, nvi, vind, yd, &
                        dy, itmax, toler, materf, expbp, &
                        taus, dalpha, dgamma, dp, crit, &
                        sgns, rp, iret)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, dgamma, msnst(1, 1, is), b_incx, gamsns, &
                       b_incy)
        end do
        do is = 1, nbsys
            call caldfp(msnst(1, 1, is), gamsns, dfpdga(1, 1, is), iret)
        end do
        nsfa = nsfa+nbsys
        nsfv = nsfv+nbsys*3
    end do
end subroutine
