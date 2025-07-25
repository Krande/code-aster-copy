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
subroutine nmisot(fami, kpg, ksp, ndim, typmod, &
                  l_epsi_varc, imate, compor, crit, deps, &
                  sigm, vim, option, sigp, vip, &
                  dsidep, iret)
!
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterf_types.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/ecpuis.h"
#include "asterfort/nmcri1.h"
#include "asterfort/nmcri2.h"
#include "asterfort/radial.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/rupmat.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/zerofr.h"
!
    aster_logical, intent(in) :: l_epsi_varc
    integer(kind=8) :: ndim, imate, kpg, ksp, iret
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: compor, option
    real(kind=8) :: crit(*), line, radi
    real(kind=8) :: deps(6), dx, deuxmu
    real(kind=8) :: sigm(6), vim(*), sigp(6), vip(*), dsidep(6, 6)
! ----------------------------------------------------------------------
!     REALISE LA LOI DE VON MISES ISOTROPE ET ELASTIQUE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!
! IN  KPG,KSP  : NUMERO DU (SOUS)POINT DE GAUSS
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
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
!----- COMMONS NECESSAIRES A VON_MISES ISOTROPE C_PLAN :
!      COMMONS COMMUNS A NMCRI1 ET NMISOT
    common/rconm1/deuxmu, nu, e, sigy, rprim, pm, sigel, line
    common/kconm1/imate2, jprol2, jvale2, nbval2
!----- COMMONS NECESSAIRES A VON_MISES ISOTROPE ECROUISSAGE PUISSANCE :
!      COMMONS COMMUNS A NMCRI2 ET NMISOT
    common/rconm2/alfafa, unsurn, sieleq
!
    aster_logical :: cplan, plasti, inco, dech
    real(kind=8) :: depsth(6), valres(3), epsthe, pm, co
    real(kind=8) :: depsmo, sigmmo, e, nu, troisk, rprim, rp, airerp
    real(kind=8) :: sieleq, sigeps, seuil, dp, coef, dsde, sigy, hydrm, hydrp
    real(kind=8) :: depsdv(6), sigmdv(6), sigpdv(6), sigdv(6)
    real(kind=8) :: em, num, troikm, deumum, sigmp(6), sigel(6), a
    real(kind=8) :: sechm, sechp, sref, tp, defam(6), defap(6)
    integer(kind=8) :: ndimsi, jprolm, jvalem, nbvalm, jprol2, jvale2, nbval2
    integer(kind=8) :: imate2, jprolp, jvalep, nbvalp, k, l, niter, ibid
    integer(kind=8) :: iret2, iret3, iret4, iret5
    integer(kind=8) :: icodre(3)
    character(len=16) :: nomres(3)
    character(len=8) :: nompar(3), para_type
    character(len=32) :: phenom
    real(kind=8) :: valpam(3), valpap(3), para_vale, valrm(2)
    real(kind=8) :: bendom, bendop, kdessm, kdessp, xm(6), xp(6)
!-----------------------------------------------------------------------
    integer(kind=8) :: lgpg
    real(kind=8) :: alfafa, coco, dp0, precr, rprim0, tm
    real(kind=8) :: unsurn, xap
!-----------------------------------------------------------------------
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    character(len=6), parameter :: epsa(6) = (/'EPSAXX', 'EPSAYY', 'EPSAZZ', 'EPSAXY', 'EPSAXZ', &
                                               'EPSAYZ'/)
! DEB ------------------------------------------------------------------
!
!     -- 1 INITIALISATIONS :
!     ----------------------
!     CES VARIABLES SONT RECOPIEES DANS LE COMMON PLUS LOIN,
!     EN LES INITIALISANT A ZERO, ON DEVRAIT VOIR ASSEZ VITE
!     SI NMCRI1 UTILISE LE COMMON SANS QUE L'ON SOIT PASSER PAR RCTRAC.
    nbvalp = 0
    jprolp = 0
    jvalep = 0
!
    cplan = typmod(1) .eq. 'C_PLAN'
    inco = typmod(2) .eq. 'INCO'
    dech = option(11:14) .eq. 'ELAS'
    if (inco) then
        co = 0.d0
    else
        co = 1.d0
    end if
    ndimsi = 2*ndim
    imate2 = imate
!
    defam = 0.d0
    defap = 0.d0
    valpam = 0.d0
    valpap = 0.d0
    valrm = 0.d0
    valres = 0.d0
    xm = 0.d0
    xp = 0.d0
    depsth = 0.d0
    depsdv = 0.d0
    sigmdv = 0.d0
    sigpdv = 0.d0
    sigdv = 0.d0
    sigmp = 0.d0
    sigel = 0.d0
    nompar = "XXXXXXXX"
!
!
    if (.not. (compor(1:9) .eq. 'VMIS_ISOT')) then
        call utmess('F', 'ALGORITH4_50', sk=compor)
    end if
!
!     -- 2 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
!    RCCOMA POUR GERER KIT_DDI (GLRC+VMIS_ISOT)
    call rccoma(imate, 'ELAS', 1, phenom, icodre(1))
    if (phenom .eq. 'ELAS') then
        nomres(1) = 'E'
        nomres(2) = 'NU'
    else if (phenom .eq. 'ELAS_GLRC') then
        nomres(1) = 'E_M'
        nomres(2) = 'NU_M'
    end if
!
    nompar(1) = 'TEMP'
    call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                ksp, tm, iret3)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, tp, iret4)
    valpam(1) = tm
    valpap(1) = tp
!
    call rcvarc(' ', 'HYDR', '-', fami, kpg, &
                ksp, hydrm, iret2)
    if (iret2 .ne. 0) hydrm = 0.d0
    call rcvarc(' ', 'HYDR', '+', fami, kpg, &
                ksp, hydrp, iret2)
    if (iret2 .ne. 0) hydrp = 0.d0
    call rcvarc(' ', 'SECH', '-', fami, kpg, &
                ksp, sechm, iret2)
    if (iret2 .ne. 0) sechm = 0.d0
    call rcvarc(' ', 'SECH', '+', fami, kpg, &
                ksp, sechp, iret2)
    if (iret2 .ne. 0) sechp = 0.d0
    call rcvarc(' ', 'SECH', 'REF', fami, kpg, &
                ksp, sref, iret2)
    if (iret2 .ne. 0) sref = 0.d0
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
    if (compor(1:14) .eq. 'VMIS_ISOT_TRAC') then
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    1, nomres(2), valres(2), icodre(2), 2)
        num = valres(2)
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    1, nomres(2), valres(2), icodre(2), 2)
        nu = valres(2)
    else
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    2, nomres(1), valres(1), icodre(1), 2)
        em = valres(1)
        num = valres(2)
        deumum = em/(1.d0+num)
!        CRIT_RUPT
        if (crit(IPOSTITER) .gt. 0.d0) then
            if (vim(8) .gt. 0.d0) then
                lgpg = 8
                call rupmat(fami, kpg, ksp, imate, vim, &
                            lgpg, em, sigm)
            end if
!       Si il y a rupture
        end if
!
        if (inco) then
            troikm = deumum
        else
            troikm = em/(1.d0-2.d0*num)
        end if
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    2, nomres(1), valres(1), icodre(1), 2)
        e = valres(1)
        nu = valres(2)
!
!        CRIT_RUPT
        if (crit(IPOSTITER) .gt. 0.d0) then
            if (vim(8) .gt. 0.d0) then
                lgpg = 8
                call rupmat(fami, kpg, ksp, imate, vim, &
                            lgpg, e, sigm)
            end if
        end if
!
        if (inco) then
            deuxmu = 2.d0*e/3.d0
            troisk = deuxmu
        else
            deuxmu = e/(1.d0+nu)
            troisk = e/(1.d0-2.d0*nu)
        end if
    end if
    call verift(fami, kpg, ksp, 'T', imate, &
                epsth_=epsthe)
!
! --- RETRAIT ENDOGENE ET RETRAIT DE DESSICCATION
!
    nomres(1) = 'B_ENDOGE'
    nomres(2) = 'K_DESSIC'
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(1), valres(1), icodre(1), 0)
    if (icodre(1) .ne. 0) valres(1) = 0.d0
    bendom = valres(1)
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(1), valres(1), icodre(1), 0)
    if (icodre(1) .ne. 0) valres(1) = 0.d0
    bendop = valres(1)
!
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(2), valres(2), icodre(2), 0)
    if (icodre(2) .ne. 0) valres(2) = 0.d0
    kdessm = valres(2)
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(2), valres(2), icodre(2), 0)
    if (icodre(2) .ne. 0) valres(2) = 0.d0
    kdessp = valres(2)
!
!     -- 3 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
    line = 0.d0
    plasti = (vim(2) .ge. 0.5d0)
    if (compor(10:14) .eq. '_LINE') then
        line = 1.d0
        nomres(1) = 'D_SIGM_EPSI'
        nomres(2) = 'SY'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 2)
        dsde = valres(1)
        sigy = valres(2)
        if ((e-dsde) .lt. r8miem()) then
            valrm(1) = dsde
            valrm(2) = e
            call utmess('F', 'COMPOR1_54', nr=2, valr=valrm)
        else
            rprim = dsde*e/(e-dsde)
        end if
        rp = rprim*vim(1)+sigy
    else if (compor(10:14) .eq. '_PUIS') then
        line = -1.d0
        nomres(1) = 'SY'
        nomres(2) = 'A_PUIS'
        nomres(3) = 'N_PUIS'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ECRO_PUIS', 0, ' ', [0.d0], &
                    3, nomres, valres, icodre, 2)
        sigy = valres(1)
        alfafa = valres(2)
        coco = e/alfafa/sigy
        unsurn = 1.d0/valres(3)
        rp = sigy*(coco*vim(1))**unsurn+sigy
        if (vim(1) .gt. r8prem()) then
            rprim = unsurn*sigy*coco*(coco*vim(1))**(unsurn-1)
        else
            rprim = e
        end if
    else if (compor(10:14) .eq. '_TRAC') then
        nompar(2) = 'SECH'
        valpam(2) = sechm
        nompar(3) = 'HYDR'
        valpam(3) = hydrm
        call rctype(imate, 3, nompar, valpam, para_vale, &
                    para_type)
!
        if ((para_type .eq. 'TEMP') .and. (iret3 .eq. 1)) then
            call utmess('F', 'COMPOR5_5', sk=para_type)
        end if
        call rctrac(imate, 1, 'SIGM', tm, jprolm, &
                    jvalem, nbvalm, em)
!
!        CRIT_RUPT VMIS_ISOT_TRAC
        if (crit(IPOSTITER) .gt. 0.d0) then
            if (vim(8) .gt. 0.d0) then
                lgpg = 8
                call rupmat(fami, kpg, ksp, imate, vim, &
                            lgpg, em, sigm)
            end if
        end if
!
        deumum = em/(1.d0+num)
        if (inco) then
            troikm = deumum
        else
            troikm = em/(1.d0-2.d0*num)
        end if
        nompar(2) = 'SECH'
        valpap(2) = sechp
        nompar(3) = 'HYDR'
        valpap(3) = hydrp
        call rctype(imate, 3, nompar, valpap, para_vale, &
                    para_type)
        if ((para_type .eq. 'TEMP') .and. (iret4 .eq. 1)) then
            call utmess('F', 'COMPOR5_5', sk=para_type)
        end if
        call rctrac(imate, 1, 'SIGM', para_vale, jprolp, &
                    jvalep, nbvalp, e)
!        CRIT_RUPT VMIS_ISOT_TRAC
        if (crit(IPOSTITER) .gt. 0.d0) then
            if (vim(8) .gt. 0.d0) then
                lgpg = 8
                call rupmat(fami, kpg, ksp, imate, vim, &
                            lgpg, e, sigm)
            end if
        end if
!
        call rcfonc('S', 1, jprolp, jvalep, nbvalp, &
                    sigy=sigy)
        call rcfonc('V', 1, jprolp, jvalep, nbvalp, &
                    p=vim(1), rp=rp, rprim=rprim, airerp=airerp)
        if (inco) then
            deuxmu = 2.d0*e/3.d0
            troisk = deuxmu
        else
            deuxmu = e/(1.d0+nu)
            troisk = e/(1.d0-2.d0*nu)
        end if
    end if
!
!     -- 4 CALCUL DE DEPSMO ET DEPSDV :
!     --------------------------------
!
    if (l_epsi_varc) then
        coef = epsthe-bendop*hydrp+bendom*hydrm-kdessp*(sref-sechp)+kdessm*(sref-sechm)
    else
        coef = 0.d0
        defam(:) = 0.d0
        defap(:) = 0.d0
    end if
    if (cplan) then
        deps(3) = -nu/(1.d0-nu)*(deps(1)+deps(2))+(1.d0+nu)/(1.d0-nu)*coef+ &
                  nu*(defap(1)-defam(1)+defap(2)-defam(2))/(1.d0-nu)+defap(3)-defam(3)
    end if
    depsmo = 0.d0
    do k = 1, 3
        depsth(k) = deps(k)-coef-(defap(k)-defam(k))
        depsmo = depsmo+depsth(k)
    end do
    depsmo = depsmo/3.d0
    do k = 4, ndimsi
        depsth(k) = deps(k)-(defap(k)-defam(k))
    end do
    do k = 1, ndimsi
        depsdv(k) = depsth(k)-depsmo*kron(k)*co
    end do
!
!     -- 5 CALCUL DE SIGMP :
!     ----------------------
    sigmmo = (sigm(1)+sigm(2)+sigm(3))/3.d0
    do k = 1, ndimsi
        sigmp(k) = deuxmu/deumum*(sigm(k)-sigmmo*kron(k))+troisk/troikm*sigmmo*kron(k)
    end do
!
!     -- 6-a CALCUL DE SIGMMO, SIGMDV, SIGEL, SIELEQ
!     --------------------------------------------
    sigmmo = (sigmp(1)+sigmp(2)+sigmp(3))/3.d0
    sieleq = 0.d0
    do k = 1, ndimsi
        sigmdv(k) = sigmp(k)-sigmmo*kron(k)
        sigel(k) = sigmdv(k)+deuxmu*depsdv(k)
        sieleq = sieleq+sigel(k)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
!
!     -- 6-b CALCUL DU SEUIL ET TRAITEMENT DU CAS DE LA QUASI-EGALITE
!     ---------------------------------------------------------------
    if (abs(sieleq-rp) .le. r8prem()*rp) then
        sieleq = rp
    end if
    seuil = sieleq-rp
!
!     -- 7 CALCUL DE SIGP,SIGPDV,VIP,DP,RP:
!     -------------------------------------
    dp = 0.d0
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
!       -- 7.1 CALCUL DE DP (ET DX SI C_PLAN) :
!       -------------------------------------------
        if (seuil .le. 0.d0) then
            vip(2) = 0.d0
            dp = 0.d0
        else
            vip(2) = 1.d0
            pm = vim(1)
            if (cplan) then
                niter = abs(nint(crit(1)))
                jprol2 = jprolp
                jvale2 = jvalep
                nbval2 = nbvalp
                precr = abs(crit(3))*sigy
!
!
!             CALCUL DE L'APPROXIMATION : DP SANS CONTRAINTE PLANE
!
                if (compor(10:14) .eq. '_LINE') then
                    dp0 = sieleq-sigy-rprim*pm
                    dp0 = dp0/(rprim+1.5d0*deuxmu)
                else if (compor(10:14) .eq. '_PUIS') then
                    dp0 = (sieleq-rp)/(1.5d0*deuxmu)
                else
                    call rcfonc('E', 1, jprolp, jvalep, nbvalp, &
                                e=e, nu=nu, p=pm, rp=rp, rprim=rprim, &
                                airerp=airerp, sieleq=sieleq, dp=dp0)
                end if
                xap = dp0
                call zerofr(0, 'DEKKER', nmcri1, 0.d0, xap, &
                            precr, niter, dp, iret, ibid)
                if (iret .eq. 1) goto 999
                if (line .gt. 0.5d0) then
                    rp = sigy+rprim*(pm+dp)
                else if (line .lt. -0.5d0) then
                    rp = sigy+sigy*(e*(pm+dp)/alfafa/sigy)**unsurn
                else
                    call rcfonc('V', 1, jprolp, jvalep, nbvalp, &
                                p=pm+dp, rp=rp, airerp=airerp)
                end if
                dx = 3.d0*(1.d0-2.d0*nu)*sigel(3)*dp/(e*dp+2.d0*(1.d0-nu)*rp)
            else
                if (compor(10:14) .eq. '_LINE') then
                    dp = sieleq-sigy-rprim*pm
                    dp = dp/(rprim+1.5d0*deuxmu)
                    rp = sigy+rprim*(pm+dp)
                else if (compor(10:14) .eq. '_PUIS') then
                    dp0 = (sieleq-rp)/(1.5d0*deuxmu)
!               AMELIORATION DE LA PREDICTION DE DP EN ESTIMANT
!               RPRIM(PM+DP0)
                    rprim0 = unsurn*sigy*coco*(coco*(pm+dp0))**(unsurn-1)
                    dp0 = dp0/(1+rprim0/1.5d0/deuxmu)
                    xap = dp0
                    precr = crit(3)*sigy
                    niter = nint(crit(1))
                    call zerofr(0, 'DEKKER', nmcri2, 0.d0, xap, &
                                precr, niter, dp, iret, ibid)
                    if (iret .eq. 1) goto 999
                    call ecpuis(e, sigy, alfafa, unsurn, pm, &
                                dp, rp, rprim)
!
                else if (compor(10:14) .eq. '_TRAC') then
                    call rcfonc('E', 1, jprolp, jvalep, nbvalp, &
                                e=e, nu=nu, p=vim(1), rp=rp, rprim=rprim, &
                                airerp=airerp, sieleq=sieleq, dp=dp)
                end if
            end if
        end if
        vip(1) = vim(1)+dp
        plasti = (vip(2) .ge. 0.5d0)
!
!
!
!
!
!         -- 7.2 CALCUL DE SIGP :
!         -----------------------
        if (cplan .and. plasti) then
            depsmo = depsmo+dx/3.d0
            depsdv(1) = depsdv(1)-dx/3.d0
            depsdv(2) = depsdv(2)-dx/3.d0
            depsdv(3) = depsdv(3)+dx*2.d0/3.d0
        end if
        do k = 1, ndimsi
            sigpdv(k) = sigmdv(k)+deuxmu*depsdv(k)
            sigpdv(k) = sigpdv(k)*rp/(rp+1.5d0*deuxmu*dp)
            sigp(k) = sigpdv(k)+(sigmmo+co*troisk*depsmo)*kron(k)
        end do
!
    end if
!
!     -- 8 CALCUL DE DSIDEP(6,6) :
!     ----------------------------
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
!
        if (option(1:10) .eq. 'RIGI_MECA_') then
!         - - OPTION='RIGI_MECA_TANG' => SIGMA(T)
            rp = 0.d0
            do k = 1, ndimsi
                sigdv(k) = sigmdv(k)
                rp = rp+sigdv(k)**2
            end do
            rp = sqrt(1.5d0*rp)
            sigeps = 1.d0
        else
!         - - OPTION='FULL_MECA' => SIGMA(T+DT)
            do k = 1, ndimsi
                sigdv(k) = sigpdv(k)
            end do
            sigeps = 0.d0
            do k = 1, ndimsi
                sigeps = sigeps+sigdv(k)*depsdv(k)
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
!      S'il YA RUPTURE ALORS INTERDIRE PLASTICITE CAR LES CONTRAINTES ONT ETE MIS A ZERO
        if (crit(IPOSTITER) .gt. 0.d0) then
            if (vim(8) .gt. 0.d0) then
                plasti = .false.
            end if
        end if
!
!
        a = 1.d0
        if (.not. dech) then
            if (plasti .and. (sigeps .ge. 0.d0)) then
                if (rp .le. 1.d-15) then
                    call utmess('F', 'ALGORITH4_46', sk=option(1:14))
                end if
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
!       -- 8.2 PARTIE ELASTIQUE:
        do k = 1, 3
            do l = 1, 3
                dsidep(k, l) = dsidep(k, l)+co*(troisk/3.d0-deuxmu/(3.d0*a))
            end do
        end do
        do k = 1, ndimsi
            dsidep(k, k) = dsidep(k, k)+deuxmu/a
        end do
!
!       -- 8.3 CORRECTION POUR LES CONTRAINTES PLANES :
        if (cplan) then
            do k = 1, ndimsi
                if (k .eq. 3) goto 136
                do l = 1, ndimsi
                    if (l .eq. 3) goto 137
                    dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep(3, l)
137                 continue
                end do
136             continue
            end do
        end if
    end if
!
    if (option(1:9) .ne. 'RIGI_MECA') then
        if (crit(10) .gt. 0.d0) then
            call radial(ndimsi, sigm, sigp, vim(2), vip(2), &
                        0, xm, xp, radi)
            if (radi .gt. crit(10)) then
                iret = 2
            end if
        end if
    end if
!
!
999 continue
!
! - "false" prediction for RIGI_MECA_TANG
!
    if (option .eq. 'RIGI_MECA_TANG') then
        sigp(1:ndimsi) = sigm(1:ndimsi)
    end if
! FIN ------------------------------------------------------------------
end subroutine
