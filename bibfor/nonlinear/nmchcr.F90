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
function nmchcr(dp)
!.======================================================================
! person_in_charge: jean-michel.proix at edf.fr
    implicit none
!
!      NMCHCR   -- CETTE ROUTINE CONCERNE L'INTEGRATION DE LA LOI
!                  DE COMPORTEMENT 'VISC_CINX_CHAB' OU 'VMIS_CINX_CHAB'
!                  RESOLUTION DE L'EQUATION SCALAIRE NON LINEAIRE EN DP
!                  (INCREMENT DE DEFORMATION PLASTIQUE CUMULEE) :
!
!  (RP/RPPMDP)*||SIGEDV- 2/3+MP1*ALPHAM1-2/3+MP1*ALPHAM2)|| = RP
!
!                  CETTE EQUATION EST RELATIVE AU MODELE DE CHABOCHE
!                  A UN OU DEUX TENSEURS CINEMATIQUES
!                  ET ELLE EST RESOLUE PAR UNE METHODE DE SECANTES
!
!   ARGUMENT        E/S  TYPE         ROLE
!    MAT(6+2*NBVAR) IN    R       TABLEAU DES COEFFICIENTS
!                                 D'ECROUISSAGE DU MATERIAU
!    DP             IN    R       INCREMENT DE DEFORMATION PLASTIQUE
!                                 CUMULEE
!    PM             IN    R       DEFORMATION PLASTIQUE CUMULEE A
!                                 L'INSTANT DU CALCUL PRECEDENT
!    NDIMSI         IN    I       DIMENSION DU VECTEUR DES CONTRAINTES
!                                 I.E. 4 EN 2D ET 6 EN 3D
!    SIGEDV(6)       IN    R       VECTEUR DES CONTRAINTES D'ESSAI, I.E.
!                                 SIGEDV = MU/(MU-)*SIGM +2MU*DELTA_EPS
!    NBVAR          IN    R       NOMBRE DE TENSEURS DE RAPPEL
!    ALFAM(6)       IN    R       LE TENSEUR DE RAPPEL XM A L'INSTANT
!    ALFA2M(6)                     DU CALCUL PRECEDENT EST RELIE
!                                 AU TENSEUR ALFAM PAR XM = 2/3*C*ALFAM
!    DEUXMU         IN    R       COEFFICIENT DE LAME :2*MU
!    VISC           IN    I       INDICATEUR DE VISCOSITE
!    MEMO           IN    I       INDICATERU EFFET MEMOIRE
!    RM             IN    R       R(INSTM)
!    RP             IN    R       R(INSTP)=RM+DR
!    QM             IN    R       Q(PM)
!    QP             OUT   R       Q(PM+DP)
!    KSIM           IN    R       KSI(PM)
!    KSIP           OUT   R       KSI(PM+DP)
!    DT             IN    R       VALEUR DE L'INCREMENT DE TEMPS DELTAT
!    F              OUT   R       VALEUR DU CRITERE DE PLASTICITE
!                                 POUR LA VALEUR DP
!
#include "asterc/r8miem.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    integer(kind=8) :: ndimsi, nbvar, visc, memo, i, idelta
    real(kind=8) :: nmchcr, dp, critme, dq, dksi(6), gq
    real(kind=8) :: epspp(6), mat(18), pm, sigedv(6), alfam(6), deuxmu
    real(kind=8) :: epspm(6), f, alfa2m(6), dt, rm, rp, qm, q, ksim(6), ksi(6)
    real(kind=8) :: r0, rinf, b, cinf, k, w, gamma0, ainf, c2inf, gamm20
    real(kind=8) :: zero, un, deux, trois, c2p, gamm2p, m2p, delta1, delta2, n1
    real(kind=8) :: n2
    real(kind=8) :: pp, cp, gammap, mp, rppmdp, seq, s(6), grjeps, norm(6)
    real(kind=8) :: mumem, valden, kvi, etam, q0mem, qmmem, dr, depsp(6)
    real(kind=8) :: rpp, coef, denom, sdenom(6), beta1, beta2
    blas_int :: b_incx, b_incy, b_n
    common/fchab/mat, pm, sigedv, epspm, alfam, alfa2m, deuxmu, rm, rp,&
     &    qm, q, ksim, ksi, dt, n1, n2, depsp,&
     &    beta1, beta2, ndimsi, nbvar, visc, memo, idelta
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     ===============
    zero = 0.0d0
    un = 1.0d0
    deux = 2.0d0
    trois = 3.0d0
!
! --- COEFFICIENTS D'ECROUISSAGE DU MATERIAU :
!     --------------------------------------
    r0 = mat(1)
    rinf = mat(2)
    b = mat(3)
    cinf = mat(4)
    k = mat(5)
    w = mat(6)
    gamma0 = mat(7)
    ainf = mat(8)
    if (nbvar .eq. 2) then
        c2inf = mat(9)
        gamm20 = mat(10)
    end if
    if (visc .eq. 1) then
        valden = mat(11)
        kvi = mat(12)
    end if
    if (memo .eq. 1) then
        etam = mat(13)
        q0mem = mat(14)
        qmmem = mat(15)
        mumem = mat(16)
    end if
    if (idelta .gt. 0) then
        delta1 = mat(17)
        delta2 = mat(18)
    else
        delta1 = 1.d0
        delta2 = 1.d0
    end if
    beta1 = 0.d0
    beta2 = 0.d0
!
! --- CALCUL DES DIFFERENTS TERMES INTERVENANT DANS LE CRITERE
! --- DE PLASTICITE :
!     =============
    pp = pm+dp
    cp = cinf*(un+(k-un)*exp(-w*pp))
    gammap = gamma0*(ainf+(un-ainf)*exp(-b*pp))
    mp = cp/(un+gammap*dp*delta1)
    if (nbvar .eq. 2) then
        c2p = c2inf*(un+(k-un)*exp(-w*pp))
        gamm2p = gamm20*(ainf+(un-ainf)*exp(-b*pp))
        m2p = c2p/(un+gamm2p*dp*delta2)
    else
        c2p = zero
        gamm2p = zero
        m2p = zero
    end if
!
! CALCUL DE LA NORMALE
    seq = zero
    do i = 1, ndimsi
        if (nbvar .eq. 1) then
            s(i) = sigedv(i)-deux/trois*mp*alfam(i)
        else if (nbvar .eq. 2) then
            s(i) = sigedv(i)-deux/trois*mp*alfam(i)-deux/trois*m2p*alfa2m(i)
        end if
        seq = seq+s(i)*s(i)
    end do
    seq = sqrt(trois/deux*seq)
    do i = 1, ndimsi
        norm(i) = sqrt(1.5d0)*s(i)/seq
    end do
!
!     R(P) SANS EFFET DE MEMOIRE
    if (memo .eq. 0) then
        rpp = rinf+(r0-rinf)*exp(-b*pp)
    end if
!
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, norm, b_incx, depsp, b_incy)
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    call dscal(b_n, dp*sqrt(1.5d0), depsp, b_incx)
!
    if (memo .eq. 1) then
!
! --- DETERMINATION DE L'INCREMENT DES DEFORMATIONS PLASTIQUES
!
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epspm, b_incx, epspp, b_incy)
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, depsp, b_incx, epspp, &
                   b_incy)
!
        grjeps = 0.0d0
        do i = 1, ndimsi
            grjeps = grjeps+(epspp(i)-ksim(i))**2
        end do
        grjeps = sqrt(grjeps*1.5d0)
        critme = grjeps/1.5d0-qm
        if (critme .le. 0.0d0) then
            dq = 0.0d0
            do i = 1, ndimsi
                dksi(i) = 0.0d0
            end do
        else
            dq = etam*critme
            coef = etam*qm+dq
            do i = 1, ndimsi
                if (coef .gt. r8miem()) then
                    dksi(i) = (1.d0-etam)*dq*(epspp(i)-ksim(i))/coef
                else
                    dksi(i) = 0.d0
                end if
            end do
!            test partie positive de <n:n*>. Utilité ?
!            NNE=0.D0
!            DO I=1,NDIMSI
!            NNE=NNE+DEPSP(I)*DKSI(I)
!            ENDDO
!            IF (NNE.LT.0.D0) THEN
!             DQ=0
!             DO i=1,NDIMSI
!             DKSI(I)=0.D0
!             KSI(I)=KSIM(I)
!             ENDDO
!            ENDIF
        end if
        q = qm+dq
        do i = 1, ndimsi
            ksi(i) = ksim(i)+dksi(i)
        end do
        gq = qmmem+(q0mem-qmmem)*exp(-2.d0*mumem*q)
        dr = b*(gq-rm)*dp/(1.d0+b*dp)
        rp = rm+dr
        rpp = r0+rp
    end if
!
!
    n1 = 1.d0
    n2 = 1.d0
    if (idelta .gt. 0) then
!        CALCUL DES BETA - N1, N2 - EFFET NON RADIAL
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        beta1 = ddot(b_n, alfam, b_incx, norm, b_incy)/sqrt(1.5d0)
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        beta2 = ddot(b_n, alfa2m, b_incx, norm, b_incy)/sqrt(1.5d0)
        if ((idelta .eq. 1) .or. (idelta .eq. 3)) then
            n1 = (1.d0+gammap*delta1*dp-gammap*(1.d0-delta1)*beta1)
            n1 = n1/(1.d0+gammap*dp)
        end if
        if ((idelta .eq. 2) .or. (idelta .eq. 3)) then
            n2 = (1.d0+gamm2p*delta2*dp-gamm2p*(1.d0-delta2)*beta2)
            n2 = n2/(1.d0+gamm2p*dp)
        end if
    end if
!
!
! POUR NORMER L'EQUATION
    denom = zero
    do i = 1, ndimsi
        if (nbvar .eq. 1) then
            sdenom(i) = sigedv(i)-deux/trois*cinf*alfam(i)
        else if (nbvar .eq. 2) then
            sdenom(i) = sigedv(i)-deux/trois*cinf*alfam(i)-deux/trois*c2inf*alfa2m(i)
        end if
        denom = denom+sdenom(i)*sdenom(i)
    end do
    denom = sqrt(trois/deux*denom)
!
    rppmdp = rpp+(trois/deux*deuxmu+mp*n1+m2p*n2)*dp
!
    if (visc .eq. 1) then
        rppmdp = rppmdp+kvi*((dp/dt)**(un/valden))
    end if
    if (denom .le. r8miem()) then
        f = seq-rppmdp
    else
        f = (seq-rppmdp)/denom
    end if
!
    nmchcr = -f
!
end function
