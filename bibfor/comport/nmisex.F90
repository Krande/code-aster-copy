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
subroutine nmisex(fami, kpg, ksp, ndim, imate, &
                  compor, crit, instam, instap, deps, &
                  sigm, vim, option, sigp, vip, &
                  typmod, dsidep, iret)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/nmcri1.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/zerofr.h"
!
    integer(kind=8) :: imate, ndim, kpg, ksp, iret
    real(kind=8) :: crit(*), instam, instap
    real(kind=8) :: deps(6), sigm(6), sigp(6), vim(*), vip(*)
    real(kind=8) :: dsidep(6, 6)
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*), option
    character(len=*) :: fami
!
! ----------------------------------------------------------------------
!     REALISE LA LOI DE VON MISES ISOTROPE ET ELASTIQUE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!     POUR LA METHODE IMPL-EX
!
! IN  KPG,KSP  : NUMERO DU (SOUS)POINT DE GAUSS
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  DEPS    : INCREMENT DE DEFORMATION
!               SI C_PLAN DEPS(3) EST EN FAIT INCONNU (ICI:0)
!                 =>  ATTENTION LA PLACE DE DEPS(3) EST ALORS UTILISEE.
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
! OUT IRET    : CODE RETOUR DE L'INTEGRATION DE LA LOI DE VOM MISES
!               = 1  => PAS DE PROBLEME
!               = 0  => ECHEC DANS L'INTEGRATION DE LA LOI
!
!
!
!
!
!
    aster_logical :: cplan, plasti, inco
    integer(kind=8) :: ndimsi, jprol2, jvale2, nbval2
    integer(kind=8) :: imate2, k, l, niter, ibid
    integer(kind=8) :: iret5
    real(kind=8) :: depsth(6), valres(3), epsthe, pm, co, dt, deuxmu
    real(kind=8) :: depsmo, sigmmo, e, nu, troisk, rprim, rp, p, dx
    real(kind=8) :: sieleq, sigeps, seuil, dp, coef, dsde, sigy
    real(kind=8) :: kron(6), depsdv(6), sigmdv(6), sigpdv(6), sigdv(6)
    real(kind=8) :: em, num, troikm, deumum, sigmp(6), sigel(6), a
    real(kind=8) :: defam(6), defap(6), line
    real(kind=8) :: precr, dp0, xap
    real(kind=8) :: valrm(2)
    real(kind=8) :: rac2
    integer(kind=8) :: icodre(3)
    character(len=6) :: epsa(6)
    character(len=16) :: nomres(3)
!
!----- COMMONS NECESSAIRES A VON_MISES ISOTROPE C_PLAN :
!      COMMONS COMMUNS A NMCRI1 ET NMISOT
    common/rconm1/deuxmu, nu, e, sigy, rprim, pm, sigel, line
    common/kconm1/imate2, jprol2, jvale2, nbval2
!
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data epsa/'EPSAXX', 'EPSAYY', 'EPSAZZ', 'EPSAXY', 'EPSAXZ',&
     &              'EPSAYZ'/
!
! DEB ------------------------------------------------------------------
!
!     -- 1 INITIALISATIONS :
!     ----------------------
    cplan = typmod(1) .eq. 'C_PLAN'
    inco = typmod(2) .eq. 'INCO'
    if (inco) then
        co = 0.d0
    else
        co = 1.d0
    end if
    dt = instap-instam
!
    ndimsi = 2*ndim
    imate2 = imate
    rac2 = sqrt(2.d0)
!
!
!
!     -- 2 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
    nomres(1) = 'E'
    nomres(2) = 'NU'
!
    defam(:) = 0.d0
    defap(:) = 0.d0
!
    do k = 1, ndimsi
        call rcvarc(' ', epsa(k), '-', fami, kpg, &
                    ksp, defam(k), iret5)
        if (iret5 .ne. 0) defam(k) = 0.d0
!
        call rcvarc(' ', epsa(k), '+', fami, kpg, &
                    ksp, defap(k), iret5)
        if (iret5 .ne. 0) defap(k) = 0.d0
    end do
!
! MISE AU FORMAT DES TERMES NON DIAGONAUX
!
    do k = 4, ndimsi
        defam(k) = defam(k)*rac2
        defap(k) = defap(k)*rac2
    end do
!
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres(1), valres(1), icodre(1), 2)
!
    em = valres(1)
    num = valres(2)
    deumum = em/(1.d0+num)
    if (inco) then
        troikm = deumum
    else
        troikm = em/(1.d0-2.d0*num)
    end if
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres(1), valres(1), icodre(1), 2)
!
    e = valres(1)
    nu = valres(2)
    deuxmu = e/(1.d0+nu)
    if (inco) then
        deuxmu = 2.d0*e/3.d0
        troisk = deuxmu
    else
        deuxmu = e/(1.d0+nu)
        troisk = e/(1.d0-2.d0*nu)
    end if
!
    call verift(fami, kpg, ksp, 'T', imate, &
                epsth_=epsthe)
!
!
!     -- 3 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
    if (compor(1) .eq. 'VMIS_ISOT_LINE') then
        plasti = (vim(2) .gt. 0.0d0)
        line = 1.d0
        nomres(1) = 'D_SIGM_EPSI'
        nomres(2) = 'SY'
!
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 2)
!
        dsde = valres(1)
        sigy = valres(2)
!
        if ((e-dsde) .lt. r8miem()) then
            valrm(1) = dsde
            valrm(2) = e
            call utmess('F', 'COMPOR1_54', nr=2, valr=valrm)
        else
            rprim = dsde*e/(e-dsde)
        end if
!
        rp = rprim*vim(1)+sigy
    else
!       -- CAS : COMPOR = 'ELAS'
        rp = 0.d0
        plasti = .false.
    end if
!
!     -- 4 CALCUL DE DEPSMO ET DEPSDV :
!     --------------------------------
    coef = epsthe
!
    if (cplan) deps(3) = -nu/(1.d0-nu)*(deps(1)+deps(2))+(1.d0+nu)/(1.d0-nu)*coef+nu*(defap(1)&
                         &-defam(1)+defap(2)-defam(2))/(1.d0-nu)+defap(3)-defam(3)
!
    depsmo = 0.d0
!
    do k = 1, 3
        depsth(k) = deps(k)-coef-(defap(k)-defam(k))
        depsth(k+3) = deps(k+3)-(defap(k+3)-defam(k+3))
        depsmo = depsmo+depsth(k)
    end do
!
    depsmo = depsmo/3.d0
!
    do k = 1, ndimsi
        depsdv(k) = depsth(k)-depsmo*kron(k)*co
    end do
!
!     -- 5 CALCUL DE SIGMP :
!     ----------------------
    sigmmo = 0.d0
!
    do k = 1, 3
        sigmmo = sigmmo+sigm(k)
    end do
!
    sigmmo = sigmmo/3.d0
!
    do k = 1, ndimsi
        sigmp(k) = deuxmu/deumum*(sigm(k)-sigmmo*kron(k))+troisk/troikm*sigmmo*kron(k)
    end do
!
!     -- 6 CALCUL DE SIGMMO, SIGMDV, SIGEL, SIELEQ ET SEUIL :
!     -------------------------------------------------------
    sigmmo = 0.d0
!
    do k = 1, 3
        sigmmo = sigmmo+sigmp(k)
    end do
!
    sigmmo = sigmmo/3.d0
    sieleq = 0.d0
!
    do k = 1, ndimsi
        sigmdv(k) = sigmp(k)-sigmmo*kron(k)
        sigel(k) = sigmdv(k)+deuxmu*depsdv(k)
        sieleq = sieleq+sigel(k)**2
    end do
!
    sieleq = sqrt(1.5d0*sieleq)
    seuil = sieleq-rp
!
!     -- 7 CALCUL DE SIGP,SIGPDV,VIP,DP,RP:
!     -------------------------------------
    dp = 0.d0
!
    if (option .eq. 'RAPH_MECA_IMPLEX') then
!
        if (compor(1) (1:4) .eq. 'ELAS') then
            do k = 1, ndimsi
                sigp(k) = sigmp(k)+deuxmu*depsdv(k)+co*troisk*depsmo*kron(k)
            end do
!
!       -- 7.1 CALCUL DE DP (ET DX SI C_PLAN) :
!       -------------------------------------------
        else if (compor(1) .eq. 'VMIS_ISOT_LINE') then
            if (seuil .le. 0.d0) then
                vip(2) = 0.d0
!            VIP(3) = 0.D0
                dp = 0.d0
            else
!            VIP(2) = 1.D0
                pm = vim(1)
!
                if (cplan) then
                    niter = abs(nint(crit(1)))
                    precr = abs(crit(3))*sigy
!
!             CALCUL DE L'APPROXIMATION : DP SANS CONTRAINTE PLANE
!
                    dp0 = sieleq-sigy-rprim*pm
                    dp0 = dp0/(rprim+1.5d0*deuxmu)
                    xap = dp0
!
                    call zerofr(0, 'DEKKER', nmcri1, 0.d0, xap, &
                                precr, niter, dp, iret, ibid)
!
                    if (iret .ne. 0) goto 999
!
                    rp = sigy+rprim*(pm+dp)
                    dx = 3.d0*(1.d0-2.d0*nu)*sigel(3)*dp/(e*dp+2.d0*(1.d0-nu)*rp)
                else
                    dp = sieleq-sigy-rprim*pm
                    dp = dp/(rprim+1.5d0*deuxmu)
                    rp = sigy+rprim*(pm+dp)
                end if
            end if
!
            vip(1) = vim(1)+dp
!          VIP(3) = DP/DT
            vip(2) = dp/dt
            plasti = (vip(2) .gt. 0.0d0)
!
!         -- 7.2 CALCUL DE SIGP :
!         -----------------------
            if (cplan .and. plasti) then
                depsmo = depsmo+dx/3.d0
                depsdv(1) = depsdv(1)-dx/3.d0
                depsdv(2) = depsdv(2)-dx/3.d0
                depsdv(3) = depsdv(3)+dx*2.d0/3.d0
            end if
!
            do k = 1, ndimsi
                sigpdv(k) = sigmdv(k)+deuxmu*depsdv(k)
                sigpdv(k) = sigpdv(k)*rp/(rp+1.5d0*deuxmu*dp)
                sigp(k) = sigpdv(k)+(sigmmo+co*troisk*depsmo)*kron(k)
            end do
        end if
!
!
!    RIGI_MECA_IMPLEX : EXTRAPOLATION, CONTRAINTES ET MATRICE
!
    else if (option .eq. 'RIGI_MECA_IMPLEX') then
!
!    EXTRAPOLATION
!
        dp = max(vim(2)*dt, 0.d0)
        p = vim(1)+dp
!
!    MISE A JOUR DE LA VARIABLE INTERNE
!
        rp = sigy+rprim*p
!
!    CONTRAINTES
!
        if (compor(1) (1:4) .eq. 'ELAS') then
            do k = 1, ndimsi
                sigp(k) = sigmp(k)+deuxmu*depsdv(k)+co*troisk*depsmo*kron(k)
            end do
!
!       -- 7.1 CALCUL DE DP (ET DX SI C_PLAN) :
!       -------------------------------------------
        else if (compor(1) .eq. 'VMIS_ISOT_LINE') then
            if (seuil .le. 0.d0) then
                dp = 0.d0
            else
                pm = vim(1)
                dp = sieleq-sigy-rprim*pm
                dp = dp/(rprim+1.5d0*deuxmu)
                rp = sigy+rprim*(pm+dp)
            end if
!
            do k = 1, ndimsi
                sigpdv(k) = sigmdv(k)+deuxmu*depsdv(k)
                sigpdv(k) = sigpdv(k)*rp/(rp+1.5d0*deuxmu*dp)
                sigp(k) = sigpdv(k)+(sigmmo+co*troisk*depsmo)*kron(k)
            end do
        end if
!
!    MATRICE TANGENTE
!
        if (option .eq. 'RIGI_MECA_IMPLEX') then
            rp = 0.d0
            do k = 1, ndimsi
                sigdv(k) = sigmdv(k)
                rp = rp+sigdv(k)**2
            end do
            rp = sqrt(1.5d0*rp)
        else
!         - - OPTION='FULL_MECA' => SIGMA(T+DT)
            if (compor(1) .eq. 'VMIS_ISOT_LINE') then
                do k = 1, ndimsi
                    sigdv(k) = sigpdv(k)
                end do
            end if
        end if
!
!       -- 8.1 PARTIE PLASTIQUE:
        do k = 1, ndimsi
            do l = 1, ndimsi
                dsidep(k, l) = 0.d0
            end do
        end do
!
        a = 1.d0
        if (compor(1) .eq. 'VMIS_ISOT_LINE') then
            sigeps = 0.d0
            do k = 1, ndimsi
                sigeps = sigeps+sigdv(k)*depsdv(k)
            end do
            if (plasti .and. sigeps .ge. 0.d0) then
                a = 1.d0+1.5d0*deuxmu*dp/rp
                coef = -(1.5d0*deuxmu)**2/(1.5d0*deuxmu+rprim)/rp**2*(1.d0-dp*rprim/rp)/a
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        dsidep(k, l) = coef*sigdv(k)*sigdv(l)
                    end do
                end do
            end if
        end if
!
!
!       -- 8.2 PARTIE ELASTIQUE:
        do k = 1, 3
            do l = 1, 3
                dsidep(k, l) = dsidep(k, l)+co*(troisk/3.d0-deuxmu/(3.d0*a))
            end do
        end do
!
        do k = 1, ndimsi
            dsidep(k, k) = dsidep(k, k)+deuxmu/a
        end do
!
!       -- 8.3 CORRECTION POUR LES CONTRAINTES PLANES :
        if (cplan) then
            do k = 1, ndimsi
                if (k .eq. 3) goto 136
!
                do l = 1, ndimsi
                    if (l .eq. 3) goto 137
!
                    dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep(3, l)
!
137                 continue
                end do
136             continue
            end do
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
999 continue
!
end subroutine
