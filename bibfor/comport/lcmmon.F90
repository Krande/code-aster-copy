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
subroutine lcmmon(fami, kpg, ksp, rela_comp, nbcomm, &
                  cpmono, nmat, nvi, vini, x, &
                  dtime, pgl, mod, coeft, neps, &
                  epsd, detot, coel, dvin, nfs, &
                  nsg, toutms, hsr, itmax, toler, &
                  iret)
! aslint: disable=W1306,W1504,W1504
    implicit none
#include "asterfort/assert.h"
#include "asterfort/calsig.h"
#include "asterfort/caltau.h"
#include "asterfort/lcgrla.h"
#include "asterfort/lcmmlc.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcrksg.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    integer(kind=8) :: kpg, ksp, nmat, nbcomm(nmat, 3), nvi, itmax, iret, nfs, nsg, neps
    real(kind=8) :: vini(*), dvin(*), x, dtime, coeft(nmat), coel(nmat)
    real(kind=8) :: sigi(6), epsd(neps), detot(neps), pgl(3, 3), toler
    real(kind=8) :: hsr(nsg, nsg, 1)
    character(len=*) :: fami
    character(len=16) :: rela_comp
!
!       IN FAMI     :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!         KPG,KSP   :  NUMERO DU (SOUS)POINT DE GAUSS
!         rela_comp   :  NOM DU MODELE DE COMPORTEMENT
!           MOD     :  TYPE DE MODELISATION
!           IMAT    :  ADRESSE DU MATERIAU CODE
!         NBCOMM :  NOMBRE DE COEF MATERIAU PAR FAMILLE
!         CPMONO :  NOMS DES LOIS MATERIAU PAR FAMILLE
!           PGL   : MATRICE DE PASSAGE GLOBAL LOCAL
!           NVI     :  NOMBRE DE VARIABLES INTERNES
!           VINI    :  VARIABLES INTERNES A T
!           X       :  INTERVALE DE TEMPS ADAPTATIF
!           DTIME   :  INTERVALE DE TEMPS
!              :  COEFFICIENTS MATERIAU INELASTIQUE A T
!           EPSD    :  DEFORMATION TOTALE A T
!           DETOT   :  INCREMENT DE DEFORMATION TOTALE
!     OUT:
!           DVIN    :  DERIVEES DES VARIABLES INTERNES A T
! INTEGRATION DES LOIS MONOCRISTALLINES PAR UNE METHODE DE RUNGE KUTTA
!
!     CETTE ROUTINE FOURNIT LA DERIVEE DE L ENSEMBLE DES VARIABLES
!     INTERNES DU MODELE
!
!     ------------------------------------------------------------------
    character(len=8) :: mod
    character(len=16) :: nomfam
    character(len=24) :: cpmono(5*nmat+1)
    real(kind=8) :: dt, dy(6+nsg), expbp(nsg), crit, sgns, q(3, 3), evi(6)
    real(kind=8) :: devi(6), mus(6), ng(3), taus, dgamma, dalpha, dp, rp
    real(kind=8) :: yd(6+nsg)
    real(kind=8) :: fkooh(6, 6), materf(nmat*2), msns(3, 3), gamsns(3, 3), lg(3)
    real(kind=8) :: toutms(nfs, nsg, 6), fp(3, 3), fp1(3, 3), deps(6), depsdt
    integer(kind=8) :: itens, nbfsys, i, nuvi, ifa, nbsys, is, nsfa, nsfv
    common/deps6/depsdt
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ------------------------------------------------------------------
! --  VARIABLES INTERNES
!
    call r8inir(9, 0.d0, gamsns, 1)
    call r8inir(nsg, 0.d0, dy, 1)
    do itens = 1, 6
        evi(itens) = vini(itens)
        devi(itens) = 0.d0
    end do
!
    b_n = to_blas_int(nmat)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, coeft, b_incx, materf(nmat+1), b_incy)
    b_n = to_blas_int(nmat)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, coel, b_incx, materf(1), b_incy)
!
!     CALCUL DU NOMBRE TOTAL DE SYSTEMES DE GLISSEMENT
    nbfsys = nbcomm(nmat, 2)
    nbsyst = 0
    do ifa = 1, nbfsys
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ng, lg, 0, q)
        nbsyst = nbsyst+nbsys
    end do
!
    if (coeft(nbcomm(1, 1)) .ge. 4) then
!         KOCKS-RAUCH ET DD_CFC : VARIABLE PRINCIPALE=DENSITE DISLOC
        ASSERT(nbcomm(nmat, 2) .eq. 1)
        do i = 1, nbsyst
            yd(6+i) = vini(6+3*(i-1)+1)
        end do
    else
!        AUTRES COMPORTEMENTS MONOCRISTALLINS
        do i = 1, nbsyst
            yd(6+i) = vini(6+3*(i-1)+2)
        end do
    end if
!
!
!     INVERSE DE L'OPERATEUR D'ELASTICITE DE HOOKE
    if (coel(nmat) .eq. 0) then
        call lcopil('ISOTROPE', mod, coel, fkooh)
    else if (coel(nmat) .eq. 1) then
        call lcopil('ORTHOTRO', mod, coel, fkooh)
    end if
!
    if (gdef .eq. 1) then
        call lcrksg(rela_comp, nvi, vini, epsd, detot, &
                    nmat, coel, sigi)
        call lcgrla(detot, deps)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, sqrt(2.d0), deps(4), b_incx)
    else
        call calsig(fami, kpg, ksp, evi, mod, &
                    rela_comp, vini, x, dtime, epsd, &
                    detot, nmat, coel, sigi)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, detot, b_incx, deps, b_incy)
    end if
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    depsdt = sqrt(ddot(b_n, deps, b_incx, deps, b_incy)/1.5d0)/dtime
    nbfsys = nbcomm(nmat, 2)
!
    nuvi = 6
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
!     NSFA : debut de la famille IFA dans DY et YD, YF
    nsfa = 6
!
    do ifa = 1, nbfsys
!
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
!
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ng, lg, 0, q)
!
        do is = 1, nbsys
!
!           CALCUL DE LA SCISSION REDUITE =
!           PROJECTION DE SIG SUR LE SYSTEME DE GLISSEMENT
!           TAU      : SCISSION REDUITE TAU=SIG:MUS
!
            call caltau(ifa, is, sigi, fkooh, nfs, &
                        nsg, toutms, taus, mus, msns)
!
!           CALCUL DE L'ECOULEMENT SUIVANT LE COMPORTEMENT
!           ECOULEMENT VISCOPLASTIQUE:
!           ROUTINE COMMUNE A L'IMPLICITE (PLASTI-LCPLNL)
!           ET L'EXPLICITE (NMVPRK-GERPAS-RK21CO-RDIF01)
!           CAS IMPLICITE : IL FAUT PRENDRE EN COMPTE DTIME
!           CAS EXPLICITE : IL NE LE FAUT PAS (ON CALCULE DES VITESSES)
!           D'OU :
            dt = -1.d0
!
            call lcmmlc(nmat, nbcomm, cpmono, nfs, nsg, &
                        hsr, nsfv, nsfa, ifa, nbsys, &
                        is, dt, nvi, vini, yd, &
                        dy, itmax, toler, materf, expbp, &
                        taus, dalpha, dgamma, dp, crit, &
                        sgns, rp, iret)
!
            if (iret .gt. 0) then
                goto 999
            end if
!
            nuvi = nuvi+3
!
            dvin(nuvi-2) = dalpha
            dvin(nuvi-1) = dgamma
            dvin(nuvi) = dp
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
! --    DERIVEES DES VARIABLES INTERNES
!
    do itens = 1, 6
        dvin(itens) = devi(itens)
    end do
    if (gdef .eq. 1) then
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vini(nvi-3-18+1), b_incx, fp, b_incy)
        fp1 = matmul(gamsns, fp)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fp1, b_incx, dvin(nvi-3-18+1), b_incy)
    end if
!
999 continue
end subroutine
