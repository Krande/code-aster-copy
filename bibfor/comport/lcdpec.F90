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
subroutine lcdpec(BEHinteg, vind, nbcomm, nmat, ndt, &
                  cpmono, materf, iter, nvi, itmax, &
                  toler, pgl, nfs, nsg, toutms, &
                  hsr, dt, dy, yd, vinf, &
                  sigf, df, nr, mod, codret)
! aslint: disable=W1306,W1504
!
    use Behaviour_type
!
    implicit none
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
!     POST-TRAITEMENTS POUR LE MONOCRISTAL
!     DEFORMATION PLASTIQUE EQUIVALENTE CUMULEE MACROSCOPIQUE
!     RECALCUL DES 3 VARIABLES INTERNES PAR SYSTEME
!  IN VIND   :  VARIABLE INTERNES A T
!     NBCOMM :  INCIDES DES COEF MATERIAU
!     NMAT   :  DIMENSION MATER ET DE NBCOMM
!     NDT    :  NOMBRE DE CMP DE SIG (6)
!     NBCOMM :  INCIDES DES COEF MATERIAU monocristal
!     MATERF :  COEF MATERIAU
!     ITER   :  NOMBRE D ITERATIONS POUR CONVERGER
!     NVI    :  NOMBRE DE VARIABLES INTERNES
!     ITMAX  :  ITER_INTE_MAXI
!     TOLER  :  RESI_INTE
!     PGL    :  MATRICE DE PASSAGE
!     TOUTMS :  TENSEURS D'ORIENTATION monocristal
!     HSR    :  MATRICE D'INTERACTION monocristal
!     DT     :  INCREMENT DE TEMPS
!     DY     :  INCREMENT DES VARIABLES Y
!     YD     :  VARIABLES A T   = ( SIGD  VARD  )
!     COMP   :  COMPOR - LOI ET TYPE DE DEFORMATION
!     SIGF   :  CONRIANTES DE CAUCHY (HPP) OU KIRCHHOFF (GDEF)
!     DF     :  GRADIENT DF
!     NR     :  DIMENSION DECLAREE DRDY
!     MOD    :  TYPE DE MODELISATION
!     CODRET :  CODE RETOUR
! VAR VINF   :  VARIABLES INTERNES A L'INSTANT ACTUEL
!
!     ----------------------------------------------------------------
#include "asterc/r8miem.h"
#include "asterfort/calcfe.h"
#include "asterfort/caltau.h"
#include "asterfort/lcgrla.h"
#include "asterfort/lcmcli.h"
#include "asterfort/lcmmlc.h"
#include "asterfort/lcmmro.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/lcnrte.h"
#include "asterfort/lcopil.h"
#include "asterfort/pk2sig.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: nmat, ndt, i, j, nbcomm(nmat, 3), nbsys, ifa, is, nbfsys, itmax
    integer(kind=8) :: nuvi, iter, nvi, iret, ir, nr, nsfa, nsfv, ifl, nuecou, codret
    integer(kind=8) :: nfs, nsg, ns, indtau, iei, is3, iv, iv3
    real(kind=8) :: vind(*), vinf(*), dy(*), materf(nmat*2)
    real(kind=8) :: epseq, pgl(3, 3), mus(6), ng(3), dgamma, dp, dalpha
    real(kind=8) :: devi(6), toutms(nfs, nsg, 6), toler, hsr(nsg, nsg)
    real(kind=8) :: taus, fkooh(6, 6), msns(3, 3), yd(*), iden(3, 3)
    real(kind=8) :: crit, sgns, dt, omp(3), qm(3, 3), fp(3, 3)
    real(kind=8) :: sicl, lg(3), rp, tau(60)
    real(kind=8) :: pk2(6), df(3, 3), id6(6), expbp(nsg)
    real(kind=8) :: fetfe6(6), gamsns(3, 3), fe(3, 3), sigf(6), rhoirr(12), xi
    real(kind=8) :: rhosat, phisat, dz, roloop(12), fivoid(12), sdp, dps(30)
    character(len=16) :: nomfam, necoul
    character(len=24) :: cpmono(5*nmat+1)
    character(len=8) :: mod
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
    data iden/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
    data id6/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
    codret = 0
    iret = 0
    sicl = -r8miem()
!     CAS MONO1 : ON RECALCULE LES VARIABLES INTERNES
    call r8inir(6, 0.d0, devi, 1)
    call r8inir(3, 0.d0, omp, 1)
!
    nbfsys = nbcomm(nmat, 2)
!
!     NSFA : debut de la famille IFA dans DY et YD
    nsfa = 6
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
!
    if (nbcomm(nmat, 1) .gt. 0) then
!        ROTATION RESEAU
        ir = 1
        do i = 1, 3
            do j = 1, 3
                qm(i, j) = vind(nvi-19+3*(i-1)+j)+iden(i, j)
            end do
        end do
    else
        ir = 0
    end if
!
    if (gdef .eq. 1) then
        if (materf(nmat) .eq. 0) then
            call lcopil('ISOTROPE', mod, materf(1), fkooh)
        else if (materf(nmat) .eq. 1) then
            call lcopil('ORTHOTRO', mod, materf(1), fkooh)
        end if
        fetfe6(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigf(1:ndt))
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        call dscal(b_n, 2.d0, fetfe6, b_incx)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, id6, b_incx, fetfe6, &
                   b_incy)
        call r8inir(9, 0.d0, gamsns, 1)
    end if
    indtau = 0
    do ifa = 1, nbfsys
!
        ifl = nbcomm(ifa, 1)
        nuecou = nint(materf(nmat+ifl))
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        necoul = cpmono(5*(ifa-1)+3) (1:16)
!
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ng, lg, 0, qm)
!
        if (necoul .eq. 'MONO_DD_CC_IRRA') then
            b_n = to_blas_int(12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vind(nsfv+3*nbsys+1), b_incx, rhoirr, b_incy)
            irr = 1
            xi = materf(nmat+ifl+23)
        else if (necoul .eq. 'MONO_DD_CFC_IRRA') then
            b_n = to_blas_int(12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vind(nsfv+3*nbsys+1), b_incx, roloop, b_incy)
            b_n = to_blas_int(12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vind(nsfv+3*nbsys+13), b_incx, fivoid, b_incy)
            irr = 2
            iei = nbcomm(ifa, 3)
            rhosat = materf(nmat+iei+8)
            phisat = materf(nmat+iei+9)
            xi = materf(nmat+iei+10)
            dz = materf(nmat+iei+11)
        else
            irr = 0
        end if
!
        do is = 1, nbsys
!
            call caltau(ifa, is, sigf, fkooh, nfs, &
                        nsg, toutms, taus, mus, msns)
!
            call lcmmlc(nmat, nbcomm, cpmono, nfs, nsg, &
                        hsr, nsfv, nsfa, ifa, nbsys, &
                        is, dt, nvi, vind, yd, &
                        dy, itmax, toler, materf, expbp, &
                        taus, dalpha, dgamma, dp, crit, &
                        sgns, rp, iret)
!
            if (iret .gt. 0) goto 999
!
            if (gdef .eq. 0) then
                do i = 1, 6
                    devi(i) = devi(i)+mus(i)*dgamma
                end do
            else
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, dgamma, msns, b_incx, gamsns, &
                           b_incy)
            end if
!
! STOCKAGE DES VARIABLES INTERNES PAR SYSTEME DE GLISSEMENT
!
            nuvi = nsfv+3*(is-1)+3
            vinf(nuvi-2) = vind(nuvi-2)+dalpha
            vinf(nuvi-1) = vind(nuvi-1)+dgamma
            vinf(nuvi) = vind(nuvi)+dp
            dps(is) = dp
            if ((nuecou .eq. 4) .or. (nuecou .eq. 5)) then
                if (vinf(nuvi-2) .lt. 0.d0) codret = 1
            end if
!
! CONTRAINTE DE CLIVAGE
            call lcmcli(nomfam, nbsys, is, pgl, sigf, &
                        sicl)
!
            call lcmmsg(nomfam, nbsys, is, pgl, mus, &
                        ng, lg, ir, qm)
            if (ir .eq. 1) then
!              ROTATION RESEAU - CALCUL DE OMEGAP
                omp(1) = omp(1)+dgamma*0.5d0*(ng(2)*lg(3)-ng(3)*lg(2))
                omp(2) = omp(2)+dgamma*0.5d0*(ng(3)*lg(1)-ng(1)*lg(3))
                omp(3) = omp(3)+dgamma*0.5d0*(ng(1)*lg(2)-ng(2)*lg(1))
            end if
!
            if (irr .eq. 1) then
                rhoirr(is) = rhoirr(is)*exp(-xi*dp)
            end if
        end do
!
        if (irr .eq. 1) then
            b_n = to_blas_int(12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rhoirr, b_incx, vinf(nsfv+3*nbsys+1), b_incy)
        end if
!
        if (irr .eq. 2) then
            do is = 1, nbsys
!              SOMME SUR COPLA(S)
                sdp = 0.d0
                do iv = 1, 12
                    is3 = (is-1)/3
                    iv3 = (iv-1)/3
!                 PARTIE POSITIVE DE ALPHA
                    if (is3 .eq. iv3) then
                        sdp = sdp+dps(iv)
                    end if
                end do
                roloop(is) = rhosat+(roloop(is)-rhosat)*exp(-xi*sdp)
                fivoid(is) = phisat+(fivoid(is)-phisat)*exp(-dz*sdp)
            end do
            b_n = to_blas_int(12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, roloop, b_incx, vinf(nsfv+3*nbsys+1), b_incy)
            b_n = to_blas_int(12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, fivoid, b_incx, vinf(nsfv+3*nbsys+13), b_incy)
        end if
        nsfa = nsfa+nbsys
        nsfv = nsfv+nbsys*3
    end do
!
    indtau = nsfv
    if (irr .eq. 1) indtau = indtau+12
    if (irr .eq. 2) indtau = indtau+24
!     CISSIONS TAU_S
    ns = 0
!     NSFA : debut de la famille IFA dans DY et YD
    nsfa = 6
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
    do ifa = 1, nbfsys
        ifl = nbcomm(ifa, 1)
        nomfam = cpmono(5*(ifa-1)+1) (1:16)
        call lcmmsg(nomfam, nbsys, 0, pgl, mus, &
                    ng, lg, 0, qm)
        do is = 1, nbsys
!           CALCUL DE LA SCISSION REDUITE =
!           PROJECTION DE SIG SUR LE SYSTEME DE GLISSEMENT
!           TAU      : SCISSION REDUITE TAU=SIG:MUS
            call caltau(ifa, is, sigf, fkooh, nfs, &
                        nsg, toutms, tau(ns+is), mus, msns)
            call lcmmlc(nmat, nbcomm, cpmono, nfs, nsg, &
                        hsr, nsfv, nsfa, ifa, nbsys, &
                        is, dt, nvi, vind, yd, &
                        dy, itmax, toler, materf, expbp, &
                        tau(ns+is), dalpha, dgamma, dp, crit, &
                        sgns, rp, iret)
        end do
        ns = ns+nbsys
        nsfa = nsfa+nbsys
        nsfv = nsfv+nbsys*3
    end do
!
    b_n = to_blas_int(ns)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, tau, b_incx, vinf(indtau+1), b_incy)
!
!     ROTATION RESEAU DEBUT
    if (ir .eq. 1) then
        call lcmmro(BEHinteg, omp, nvi, vind, vinf)
    end if
! ROTATION RESEAU FIN
!
    if (gdef .eq. 1) then
        call calcfe(nr, ndt, nvi, vind, df, &
                    gamsns, fe, fp, iret)
        if (iret .gt. 0) goto 999
!
!        CALCUL DES CONTRAINTES DE KIRCHOFF
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigf, b_incx, pk2, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, sqrt(2.d0), pk2(4), b_incx)
        call pk2sig(3, fe, 1.d0, pk2, sigf, &
                    1)
!
! les racine(2) attendues par NMCOMP :-)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, sqrt(2.d0), sigf(4), b_incx)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, iden, b_incx, fe, &
                   b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fe, b_incx, vinf(nvi-3-18+10), b_incy)
!
        call lcgrla(fp, devi)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, devi, b_incx, vinf, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, sqrt(2.d0), devi(4), b_incx)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, iden, b_incx, fp, &
                   b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fp, b_incx, vinf(nvi-3-18+1), b_incy)
!
        epseq = lcnrte(devi)
        vinf(nvi-1) = epseq
!
    else
        do i = 1, 6
            vinf(i) = vind(i)+devi(i)
        end do
        epseq = lcnrte(devi)
        vinf(nvi-1) = vind(nvi-1)+epseq
    end if
!
    vinf(nvi-2) = sicl
!
    vinf(nvi) = iter
!
999 continue
    codret = max(codret, iret)
end subroutine
