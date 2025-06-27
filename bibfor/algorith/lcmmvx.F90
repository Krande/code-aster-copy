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
subroutine lcmmvx(sigf, vin, nmat, materf, nbcomm, &
                  cpmono, pgl, nvi, hsr, nfs, &
                  nsg, toutms, timed, timef, deps, &
                  seuil)
! aslint: disable=W1306
    implicit none
!     MONOCRISTAL  :  CALCUL DU SEUIL POUR MONOCRISTAL
!     ----------------------------------------------------------------
!     IN  FAMI   :  FAMILLE DE POINTS DE GAUSS
!     IN  KPG    :  NUMERO DU POINT DE GAUSS
!     IN  KSP    :  NUMERO DU SOUS-POINT DE GAUSS
!     IN  SIGF   :  CONTRAINTE
!     IN  VIN    :  VARIABLES INTERNES = ( X1 X2 P )
!     IN  NMAT   :  DIMENSION MATER
!     IN  MATERF :  COEFFICIENTS MATERIAU A TEMP
!         COMP   :  NOM COMPORTEMENT
!         NBCOMM :  INCIDES DES COEF MATERIAU
!         CPMONO :  NOM DES COMPORTEMENTS
!         PGL    :  MATRICE DE PASSAGE
!         NR     :  DIMENSION DECLAREE DRDY
!         NVI    :  NOMBRE DE VARIABLES INTERNES
!         HSR    :  MATRICE D'INTERACTION
!         TOUTMS :  TENSEURS D'ORIENTATION
!     OUT SEUIL  :  SEUIL  ELASTICITE
!     ----------------------------------------------------------------
#include "asterfort/lcmmfe.h"
#include "asterfort/lcmmfi.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
    integer(kind=8) :: nmat, nvi, nsfa, nsfv, iexp, nfs, nsg
    integer(kind=8) :: nbfsys, i, nuvi, ifa, nbsys, is
    integer(kind=8) :: nbcomm(nmat, 3), iret
    real(kind=8) :: sigf(6), vin(nvi), rp, hsr(nsg, nsg), deps(6)
    real(kind=8) :: materf(nmat*2), seuil, dt, dy(nvi), alpham
    real(kind=8) :: ms(6), ng(3), q(3, 3), timed, timef, lg(3), depsdt
    real(kind=8) :: taus, dgamma, dalpha, dp, expbp(nsg), depst(6)
    real(kind=8) :: pgl(3, 3), crit, sgns, toutms(nfs, nsg, 6), gammam
    character(len=24) :: cpmono(5*nmat+1)
    character(len=16) :: nomfam, necoul, necris
    common/deps6/depsdt
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
!
    seuil = -1.d0
    dt = timef-timed
    nbfsys = nbcomm(nmat, 2)
    call r8inir(nvi, 0.d0, dy, 1)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, deps, b_incx, depst, b_incy)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    depsdt = sqrt(ddot(b_n, depst, b_incx, depst, b_incy)/1.5d0)/dt
!
!     NSFV : debut de la famille IFA dans les variables internes
    nsfv = 6
!     NSFA : debut de la famille IFA dans DY et YD, YF
    nsfa = 6
    do ifa = 1, nbfsys
!
        nomfam = cpmono(5*(ifa-1)+1)
        necoul = cpmono(5*(ifa-1)+3)
        necris = cpmono(5*(ifa-1)+4)
!
        call lcmmsg(nomfam, nbsys, 0, pgl, ms, &
                    ng, lg, 0, q)
!
        if (nbsys .eq. 0) then
            call utmess('F', 'ALGORITH_70')
        end if
!
        do is = 1, nbsys
!
            nuvi = nsfv+3*(is-1)
            alpham = vin(nuvi+1)
            gammam = vin(nuvi+2)
!
!           CALCUL DE LA SCISSION REDUITE =
!           PROJECTION DE SIG SUR LE SYSTEME DE GLISSEMENT
!           TAU      : SCISSION REDUITE TAU=SIG:MS
            do i = 1, 6
                ms(i) = toutms(ifa, is, i)
            end do
!
            taus = 0.d0
            do i = 1, 6
                taus = taus+sigf(i)*ms(i)
            end do
!
!           ECROUISSAGE ISOTROPE
!
            if (necoul .ne. 'MONO_DD_KR') then
                iexp = 0
                if (is .eq. 1) iexp = 1
                call lcmmfi(materf(nmat+1), ifa, nmat, nbcomm, necris, &
                            is, nbsys, vin, nsfv, dy(nsfa+1), &
                            nfs, nsg, hsr, iexp, expbp, &
                            rp)
            end if
!
!           ECOULEMENT VISCOPLASTIQUE
!
            decal = nsfv
            call lcmmfe(taus, materf(nmat+1), materf, ifa, nmat, &
                        nbcomm, necoul, is, nbsys, vin, &
                        dy(nsfa+1), rp, alpham, gammam, dt, &
                        dalpha, dgamma, dp, crit, sgns, &
                        nfs, nsg, hsr, iret)
!
            if (iret .gt. 0) then
                dp = 1.d0
            end if
            if (dp .gt. 0.d0) then
                seuil = 1.d0
                goto 999
            end if
!
        end do
!
        nsfa = nsfa+nbsys
        nsfv = nsfv+3*nbsys
    end do
999 continue
end subroutine
