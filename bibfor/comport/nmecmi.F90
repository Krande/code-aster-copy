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
subroutine nmecmi(fami, kpg, ksp, ndim, typmod, &
                  imate, compor, crit, deps, sigm, &
                  vim, option, sigp, vip, dsidep, &
                  iret)
!
! aslint: disable=
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmcri5.h"
#include "asterfort/radial.h"
#include "asterfort/rcfon2.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/zerofr.h"
!
    integer(kind=8) :: kpg, ksp, ndim, imate, iret, iret1, iret2
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*), option
    real(kind=8) :: crit(10), tp2, line, radi
    real(kind=8) :: deps(6), prec, dx, deuxmu
    real(kind=8) :: sigm(6), vim(8), sigp(6), vip(8), dsidep(6, 6)
! ----------------------------------------------------------------------
!     REALISE LA LOI DE VON MISES ISOTROPE ET ELASTIQUE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!
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
! OUT IRET    : CODE RETOUR DE  L'INTEGRATION DE LA LOI DE VOM MISES
!                   IRET=0 => PAS DE PROBLEME
!                   IRET=1 => ABSENCE DE CONVERGENCE LORS DE
!                                        L'INTEGRATION DE LA LOI
!
!----- COMMONS NECESSAIRES A VON_MISES ISOTROPE C_PLAN :
!      COMMONS COMMUNS A NMCRI1 ET NMECMI
    common/rconm5/deuxmu, troisk, sigy, rprim, pm, sigel, tp2, line, prag, xm
    common/kconm1/imate2, jprol2, jvale2, nbval2
!
!
!
!
    aster_logical :: cplan, plasti
    real(kind=8) :: depsth(6), valres(3), pm, xp(6), plast, para_vale
    real(kind=8) :: depsmo, sigmmo, e, nu, troisk, rprim, rp, hp, gp, g1, rpm
    real(kind=8) :: sieleq, sigeps, seuil, dp, coef, dsde, sigy, xm(6)
    real(kind=8) :: sigedv(6)
    real(kind=8) :: depsdv(6), sigmdv(6), sigpdv(6), sigdv(6), cc
    real(kind=8) :: em, num, troikm, deumum, sigmp(6), sigel(6)
    real(kind=8) :: hsg, pp, prag, pragm, precr, tm, tp, epsthe
    integer(kind=8) :: ndimsi, jprolm, jvalem, nbvalm, jprol2, jvale2, nbval2
    integer(kind=8) :: jprolp, jvalep, nbvalp, k, l, niter, imate2, ibid
    integer(kind=8) :: icodre(3)
    character(len=16) :: nomres(3)
    character(len=8) :: para_type
!-----------------------------------------------------------------------
    real(kind=8) :: dp0, xap
!-----------------------------------------------------------------------
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
! DEB ------------------------------------------------------------------
!
!     -- 1 INITIALISATIONS :
!     ----------------------
    cplan = typmod(1) .eq. 'C_PLAN'
    ndimsi = 2*ndim
    imate2 = imate
    iret = 0
    jprolp = 1
!
! MISE AU FORMAT DES CONTRAINTES DE RAPPEL
!
    pm = vim(1)
    plast = vim(2)
    do k = 1, 3
        xm(k) = vim(k+2)
    end do
    do k = 4, ndimsi
        xm(k) = vim(k+2)*sqrt(2.d0)
    end do
    dp = 0.d0
!
!
    if (.not. (compor(1) (1:9) .eq. 'VMIS_ECMI')) then
        call utmess('F', 'ALGORITH4_50', sk=compor(1))
    end if
!
!
!
!     -- 2 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
!
    if (compor(1) (1:14) .eq. 'VMIS_ECMI_TRAC') then
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, nomres(2), valres(2), icodre(2), 2)
        num = valres(2)
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, nomres(2), valres(2), icodre(2), 2)
        nu = valres(2)
    else
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres(1), valres(1), icodre(1), 2)
        em = valres(1)
        num = valres(2)
        deumum = em/(1.d0+num)
        troikm = em/(1.d0-2.d0*num)
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres(1), valres(1), icodre(1), 2)
        e = valres(1)
        nu = valres(2)
        deuxmu = e/(1.d0+nu)
        troisk = e/(1.d0-2.d0*nu)
    end if
    call verift(fami, kpg, ksp, 'T', imate, &
                epsth_=epsthe)
!
!     -- 3 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
    nomres(1) = 'C'
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'PRAGER', 0, ' ', [0.d0], &
                1, nomres, valres, icodre, 1)
    prag = valres(1)
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'PRAGER', 0, ' ', [0.d0], &
                1, nomres, valres, icodre, 1)
    pragm = valres(1)
    line = 0.d0
    if (compor(1) (10:14) .eq. '_LINE') then
        line = 1.d0
        nomres(1) = 'D_SIGM_EPSI'
        nomres(2) = 'SY'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 1)
        dsde = valres(1)
        sigy = valres(2)
        rprim = dsde*e/(e-dsde)-1.5d0*prag
        rpm = rprim*pm+sigy
    else
!
        call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                    ksp, tm, iret2)
        call rctype(imate, 1, 'TEMP', [tm], para_vale, &
                    para_type)
        if ((para_type .eq. 'TEMP') .and. (iret2 .eq. 1)) then
            call utmess('F', 'COMPOR5_5', sk=para_type)
        end if
        call rctrac(imate, 1, 'SIGM', para_vale, jprolm, &
                    jvalem, nbvalm, em)
        deumum = em/(1.d0+num)
        troikm = em/(1.d0-2.d0*num)
!
        call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                    ksp, tp, iret1)
        call rctype(imate, 1, 'TEMP', [tp], para_vale, &
                    para_type)
        if ((para_type .eq. 'TEMP') .and. (iret1 .eq. 1)) then
            call utmess('F', 'COMPOR5_5', sk=para_type)
        end if
        call rctrac(imate, 1, 'SIGM', para_vale, jprolp, &
                    jvalep, nbvalp, e)
        deuxmu = e/(1.d0+nu)
        troisk = e/(1.d0-2.d0*nu)
!
        call rcfon2('S', jprolp, jvalep, nbvalp, sigy=sigy)
        call rcfon2('V', jprolp, jvalep, nbvalp, p=pm, &
                    rp=rpm, rprim=rprim, c=prag)
    end if
!
!     -- 4 CALCUL DE DEPSMO ET DEPSDV :
!     --------------------------------
    coef = epsthe
    if (cplan) deps(3) = -nu/(1.d0-nu)*(deps(1)+deps(2))+(1.d0+nu)/(1.d0-nu)*coef
    depsmo = 0.d0
    do k = 1, 3
        depsth(k) = deps(k)-coef
        depsmo = depsmo+depsth(k)
    end do
    depsmo = depsmo/3.d0
    do k = 4, ndimsi
        depsth(k) = deps(k)
    end do
    do k = 1, ndimsi
        depsdv(k) = depsth(k)-depsmo*kron(k)
    end do
!
!     -- 5 CALCUL DE SIGMP :
!     ----------------------
    sigmmo = 0.d0
    do k = 1, 3
        sigmmo = sigmmo+sigm(k)
    end do
    sigmmo = sigmmo/3.d0
    do k = 1, ndimsi
        sigmp(k) = deuxmu/deumum*(sigm(k)-sigmmo*kron(k))+troisk/ &
                   troikm*sigmmo*kron(k)
    end do
!
!     -- 6 CALCUL DE SIGMMO, SIGMDV, SIGEL, SIELEQ ET SEUIL :
!     -------------------------------------------------------
    sigmmo = 0.d0
    do k = 1, 3
        sigmmo = sigmmo+sigmp(k)
    end do
    sigmmo = sigmmo/3.d0
    sieleq = 0.d0
    do k = 1, ndimsi
        sigmdv(k) = sigmp(k)-sigmmo*kron(k)
        if (pragm .ne. 0.d0) then
            xm(k) = xm(k)*prag/pragm
        end if
        sigel(k) = sigmdv(k)+deuxmu*depsdv(k)-xm(k)
        sieleq = sieleq+sigel(k)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
    seuil = sieleq-rpm
    hp = 1.d0
    gp = 1.d0
!
!     -- 7 CALCUL DE SIGP,SIGPDV,VIP,DP,RP:
!     -------------------------------------
!
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
!       -- 7.1 CALCUL DE DP (ET DX SI C_PLAN) :
!       -------------------------------------------
        if (seuil .le. 0.d0) then
            plast = 0.d0
            dp = 0.d0
            rp = rpm
        else
            plast = 1.d0
            if (cplan) then
                prec = abs(crit(3))
                niter = abs(nint(crit(1)))
                precr = prec*sigy
!
!             CALCUL DE L'APPROXIMATION : DP SANS CONTRAINTE PLANE
!
                if (compor(1) (10:14) .eq. '_LINE') then
                    dp0 = sieleq-sigy-rprim*pm
                    dp0 = dp0/(rprim+1.5d0*(deuxmu+prag))
                else
                    jprol2 = jprolp
                    jvale2 = jvalep
                    nbval2 = nbvalp
                    call rcfon2('E', jprolp, jvalep, nbvalp, e=e, &
                                nu=nu, p=pm, rp=rp, rprim=rprim, c=prag, &
                                sieleq=sieleq, dp=dp0)
                end if
                xap = dp0
                call zerofr(0, 'DEKKER', nmcri5, 0.d0, xap, &
                            precr, niter, dp, iret, ibid)
                if (iret .eq. 1) goto 999
                if (line .ge. 0.5d0) then
                    rp = sigy+rprim*(pm+dp)
                else
                    call rcfon2('V', jprolp, jvalep, nbvalp, p=pm+dp, &
                                rp=rp, rprim=rprim, c=prag)
                end if
            else
                if (compor(1) (10:14) .eq. '_LINE') then
                    dp = sieleq-sigy-rprim*pm
                    dp = dp/(rprim+1.5d0*(deuxmu+prag))
                    rp = sigy+rprim*(pm+dp)
                else
                    call rcfon2('E', jprolp, jvalep, nbvalp, e=e, &
                                nu=nu, p=pm, rp=rp, rprim=rprim, c=prag, &
                                sieleq=sieleq, dp=dp)
                end if
            end if
        end if
        pp = pm+dp
        gp = 1.d0+1.5d0*prag*dp/rp
        hp = gp+1.5d0*deuxmu*dp/rp
        plasti = (plast .ge. 0.5d0)
!
!         -- 7.2 CALCUL DE SIGP :
!         -----------------------
        if (cplan .and. plasti) then
            hsg = hp/gp
            dx = (hsg-1.d0)*sigel(3)
            dx = dx/(deuxmu/1.5d0+troisk*hsg/3.d0)
            depsmo = depsmo+dx/3.d0
            depsdv(1) = depsdv(1)-dx/3.d0
            depsdv(2) = depsdv(2)-dx/3.d0
            depsdv(3) = depsdv(3)+dx*2.d0/3.d0
        end if
        do k = 1, ndimsi
            sigedv(k) = sigmdv(k)+deuxmu*depsdv(k)
            g1 = 1.5d0*prag*dp/rp/hp
            xp(k) = xm(k)*(1.d0-g1)+g1*sigedv(k)
            sigpdv(k) = sigedv(k)*gp/hp+xm(k)*1.5d0*deuxmu*dp/rp/hp
            sigp(k) = sigpdv(k)+(sigmmo+troisk*depsmo)*kron(k)
        end do
    end if
!
!     -- 8 CALCUL DE DSIDEP(6,6) :
!     ----------------------------
    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
!
        plasti = (plast .ge. 0.5d0)
        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
!         - - OPTION='RIGI_MECA_TANG' => SIGMA(T)
            do k = 1, ndimsi
                sigdv(k) = sigmdv(k)-xm(k)
            end do
            rp = rpm
        else
!         - - OPTION='FULL_MECA' => SIGMA(T+DT)
            do k = 1, ndimsi
                sigdv(k) = sigpdv(k)-xp(k)
            end do
        end if
!
!       -- 8.1 PARTIE PLASTIQUE:
        do k = 1, ndimsi
            do l = 1, ndimsi
                dsidep(k, l) = 0.d0
            end do
        end do
!
        sigeps = 0.d0
        do k = 1, ndimsi
            sigeps = sigeps+sigdv(k)*depsdv(k)
        end do
        if (plasti .and. sigeps .ge. 0.d0) then
            cc = -(1.5d0*deuxmu)**2/(1.5d0*(deuxmu+prag)+rprim)/rp**2*(1.d0-dp*rprim/rp)/hp
            do k = 1, ndimsi
                do l = 1, ndimsi
                    dsidep(k, l) = cc*sigdv(k)*sigdv(l)
                end do
            end do
        end if
!
!       -- 8.2 PARTIE ELASTIQUE:
        do k = 1, 3
            do l = 1, 3
                dsidep(k, l) = dsidep(k, l)+troisk/3.d0-deuxmu/3.d0*gp/hp
            end do
        end do
        do k = 1, ndimsi
            dsidep(k, k) = dsidep(k, k)+deuxmu*gp/hp
        end do
!
!       -- 8.3 CORRECTION POUR LES CONTRAINTES PLANES :
        if (cplan) then
            do k = 1, ndimsi
                if (k .eq. 3) goto 136
                do l = 1, ndimsi
                    if (l .eq. 3) goto 137
                    dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep( &
                                   k, 3)*dsidep(3, l)
137                 continue
                end do
136             continue
            end do
        end if
    end if
!
    if (option(1:9) .ne. 'RIGI_MECA') then
        if (crit(10) .gt. 0.d0) then
            call radial(ndimsi, sigm, sigp, vim(2), plast, &
                        1, vim(3), vip(3), radi)
            if (radi .gt. crit(10)) then
                iret = 2
            end if
        end if
    end if
!
! MISE AU FORMAT DES CONTRAINTES DE RAPPEL
!
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        vip(1) = pp
        vip(2) = plast
        do k = 1, 3
            vip(k+2) = xp(k)
        end do
        do k = 4, ndimsi
            vip(k+2) = xp(k)/sqrt(2.d0)
        end do
    end if
!
999 continue
end subroutine
