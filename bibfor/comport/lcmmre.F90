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
subroutine lcmmre(typmod, nmat, materd, materf, nbcomm, &
                  cpmono, pgl, nfs, nsg, toutms, &
                  hsr, nr, nvi, vind, itmax, &
                  toler, timed, timef, yd, yf, &
                  deps, dy, r, iret)
! aslint: disable=W1306,W1504
    implicit none
!     MONOCRISTAL  : CALCUL DES RESIDUS DU SYSTEME NL A RESOUDRE = R(DY)
!                    CF. R5.03.11
!                    DY =  DSIG (DBETA PAR SYSTEME)
!                    Y  =  SIG   (BETA  PAR SYSTEME)
!                    R  = ( R1  R2   )
!                    ATTENTION IL FAUT CALCULER -R
!
!     IN  TYPMOD :  TYPE DE MODELISATION
!         NMAT   :  DIMENSION MATER
!         MATERD :  COEFFICIENTS MATERIAU A T
!         MATERF :  COEFFICIENTS MATERIAU A T+DT
!         COMP   :  NOM COMPORTEMENT
!         NBCOMM :  INCIDES DES COEF MATERIAU
!         CPMONO :  NOM DES COMPORTEMENTS
!         PGL    :  MATRICE DE PASSAGE
!         TOUTMS :  TENSEURS D'ORIENTATION
!         HSR    :  MATRICE D'INTERACTION
!         NR     :  DIMENSION DECLAREE DRDY
!         NVI    :  NOMBRE DE VARIABLES INTERNES
!         VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!         ITMAX  :  ITER_INTE_MAXI
!         TOLER  :  RESI_INTE
!         TIMED  :  ISTANT PRECEDENT
!         TIMEF  :  INSTANT ACTUEL
!         YD     :  VARIABLES A T       = ( SIGD BETAD )
!         YF     :  VARIABLES A T + DT  = ( SIGF BETAF )
!         DEPS   :  INCREMENT DE DEFORMATION OU GRADIENT DF
!         DY     :  SOLUTION  =  ( DSIG DBETA )
!         NR     :  DIMENSION DECLAREE DRDY
!     OUT R      :  RESIDU DU SYSTEME NL A T + DT
!         IRET   :  CODE RETOUR
!     ----------------------------------------------------------------
#include "asterfort/calcfe.h"
#include "asterfort/caltau.h"
#include "asterfort/lcgrla.h"
#include "asterfort/lcmmlc.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/lcopil.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    integer(kind=8) :: ndt, ndi, nmat, nr, nvi, nsfv, iret
    integer(kind=8) :: nbfsys, ifa, nbsys, is, itmax, nfs, nsg
    integer(kind=8) :: nbcomm(nmat, 3), nsfa, ifl, nuecou
!
    real(kind=8) :: dkooh(6, 6), fkooh(6, 6), timed, timef
    real(kind=8) :: sigf(6), sigd(6), msns(3, 3), pgl(3, 3), dgamm1
    real(kind=8) :: deps(*), depse(6), devi(6), dt
    real(kind=8) :: epsed(6), epsgl(6), h1sigf(6), vind(*)
    real(kind=8) :: materd(nmat*2), materf(nmat*2), epsef(6)
    real(kind=8) :: mus(6), ng(3), taus, dgamma, dalpha, dp, rp, depsdt
    real(kind=8) :: r(nr), dy(nr), yd(nr), yf(nr), toler, fe(3, 3)
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg), q(3, 3), lg(3)
    real(kind=8) :: gamsns(3, 3), fp(3, 3), depst(6)
    real(kind=8) :: crit, sgns, expbp(nsg)
    character(len=8) :: typmod
    character(len=16) :: nomfam
    character(len=24) :: cpmono(5*nmat+1)
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/deps6/depsdt
!     ----------------------------------------------------------------
!
    dt = timef-timed
!     INVERSE DE L'OPERATEUR D'ELASTICITE DE HOOKE
    if (materf(nmat) .eq. 0) then
        call lcopil('ISOTROPE', typmod, materd(1), dkooh)
        call lcopil('ISOTROPE', typmod, materf(1), fkooh)
    else if (materf(nmat) .eq. 1) then
        call lcopil('ORTHOTRO', typmod, materd(1), dkooh)
        call lcopil('ORTHOTRO', typmod, materf(1), fkooh)
    end if
!
    call r8inir(9, 0.d0, gamsns, 1)
    sigf(1:ndt) = yf(1:ndt)
    call r8inir(6, 0.d0, devi, 1)
!
!     POUR DD_CC
    if (gdef .eq. 1) then
        call lcgrla(deps, depst)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, sqrt(2.d0), depst(4), b_incx)
    else
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deps, b_incx, depst, b_incy)
    end if
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    depsdt = sqrt(ddot(b_n, depst, b_incx, depst, b_incy)/1.5d0)/dt
!
    nbfsys = nbcomm(nmat, 2)
    iret = 0
!
!
!     NSFA : debut de la famille IFA dans DY et YD, YF
    nsfa = 6
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
!
    do ifa = 1, nbfsys
!
        ifl = nbcomm(ifa, 1)
        nuecou = nint(materf(nmat+ifl))
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
!
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ng, lg, 0, q)
!
        do is = 1, nbsys
!           CALCUL DE LA SCISSION REDUITE
            call caltau(ifa, is, sigf, fkooh, nfs, &
                        nsg, toutms, taus, mus, msns)
!           CALCUL DE L'ECOULEMENT SUIVANT LE COMPORTEMENT
            call lcmmlc(nmat, nbcomm, cpmono, nfs, nsg, &
                        hsr, nsfv, nsfa, ifa, nbsys, &
                        is, dt, nvi, vind, yd, &
                        dy, itmax, toler, materf, expbp, &
                        taus, dalpha, dgamma, dp, crit, &
                        sgns, rp, iret)
!
            if (iret .gt. 0) then
                goto 999
            end if
!
            if (nuecou .ge. 4) then
!           POUR LES LOIS DD_* ALPHA représente la variable principale
                r(nsfa+is) = -(dy(nsfa+is)-dalpha)
            else
                dgamm1 = dy(nsfa+is)
                r(nsfa+is) = -(dgamm1-dgamma)
            end if
!
            if (gdef .eq. 0) then
                b_n = to_blas_int(6)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, dgamma, mus, b_incx, devi, &
                           b_incy)
            else
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, dgamma, msns, b_incx, gamsns, &
                           b_incy)
            end if
        end do
!
        nsfa = nsfa+nbsys
        nsfv = nsfv+nbsys*3
!
    end do
!
    if (gdef .eq. 1) then
        call calcfe(nr, ndt, nvi, vind, deps, &
                    gamsns, fe, fp, iret)
        if (iret .gt. 0) then
            goto 999
        end if
        call lcgrla(fe, epsgl)
        h1sigf(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigf(1:ndt))
        r(1:ndt) = epsgl(1:ndt)-h1sigf(1:ndt)
    else
        sigd(1:ndt) = yd(1:ndt)
        epsed(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigd(1:ndt))
        depse(1:ndt) = deps(1:ndt)-devi(1:ndt)
        epsef(1:ndt) = epsed(1:ndt)+depse(1:ndt)
! LA PREMIERE EQUATION EST  (HF-1)SIGF -(HD-1)SIGD -(DEPS-DEPSP)=0
        h1sigf(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigf(1:ndt))
        r(1:ndt) = epsef(1:ndt)-h1sigf(1:ndt)
    end if
!
999 continue
end subroutine
