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
! aslint: disable=W1501
!
subroutine nmccam(fami, kpg, ksp, ndim, &
                  typmod, imate, carcri, &
                  deps, sigm, pcrm, option, sigp, &
                  pcrp, dsidep, retcom)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/mgauss.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/get_varc.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8) :: ndim, imate, retcom
    character(len=8) :: typmod(*)
    character(len=16) :: option
    real(kind=8) :: carcri(3)
    real(kind=8) :: deps(6), deuxmu
    real(kind=8) :: sigm(6), pcrm(7), sigp(6), pcrp(7), dsidep(6, 6)
! ----------------------------------------------------------------------
!     REALISE LA LOI DE CAM CLAY ELASTOPLASTIQUE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  carcri  : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  TM      : TEMPERATURE A L'INSTANT PRECEDENT
! IN  TP      : TEMPERATURE A L'INSTANT DU CALCUL
! IN  TREF    : TEMPERATURE DE REFERENCE
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  PCRM    : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT PCRP    : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!
!
!
!
    aster_logical :: cplan, lVari
    integer(kind=8) :: iadzi, iazk24, iret
    real(kind=8), parameter :: epxmax = 5.d0
    character(len=8) :: nomail
    real(kind=8) :: valres(10), valpam(3)
    real(kind=8) :: mu, lambda, kapa, poro, prescr, m, kcam, ptrac
    real(kind=8) :: coef, young, nu, e0, xk0, xk, fonc
    real(kind=8) :: depsmo, deppmo, depseq
    real(kind=8) :: sigmmo, sigpmo, deltap, sieqm, sieqp, sieleq, simoel, spards
    real(kind=8) :: kron(6), depsdv(6), depsth(6), sigmdv(6), sigpdv(6)
    real(kind=8) :: deltas(6), sigel(6), tplus(6)
    real(kind=8) :: a(6), aa(6), fv(6)
    real(kind=8) :: ffi(6, 6), ee(6, 6), c(6, 6), cc(6, 6)
    real(kind=8) :: v(6, 6), s(6, 6), t(6, 6), vv(6, 6)
    real(kind=8) :: hh(6, 6), ses(6, 6), gg(6, 6), sps(6, 6), hhm(6, 6)
    real(kind=8) :: d1g(6, 6), d1ghhm(6, 6), id2(6, 6), devhyd(6, 6)
    real(kind=8) :: devhym(6, 6)
    real(kind=8) :: f1, f2, f3, f4, f, fp
    real(kind=8) :: f1p, f2p, f3p, f4p
    real(kind=8) :: fxi1, fxi2, fxi3, fxi4, fxi
    real(kind=8) :: h, hp, xc, xd, xlam, xa, xu, xg, xh, xe, xf, xv, xi, rap
    real(kind=8) :: ct, v0, seuil
    real(kind=8) :: xinf, xsup, rbid
    real(kind=8) :: diff2
    real(kind=8) :: zero, un, deux, trois, six, unsde, tol, ptit
    real(kind=8) :: valm, valp, tm, tp, tref, pcrpLoc(7)
    integer(kind=8) :: ndimsi, signf, signfi
    integer(kind=8) :: i, k, l, iter, matr
    integer(kind=8) :: icodre(9)
    character(len=18) :: nomres(10)
    character(len=8) :: nompar(10)
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    parameter(six=6.d0)
    parameter(unsde=0.5d0)
!
!
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data tol/1.d-6/
! DEB ------------------------------------------------------------------
!
!     -- 1 INITIALISATIONS :
!     ----------------------
    lVari = L_VARI(option)
    cplan = typmod(1) .eq. 'C_PLAN'
    ndimsi = 2*ndim
    retcom = 0
    sigp = 0.d0
    dsidep = 0.d0
!
    ptit = r8miem()
!
! - Get temperatures
!
    call get_varc(fami, kpg, ksp, 'T', &
                  tm, tp, tref)
!
!     -- 2 RECUPERATION DES CARACTERISTIQUES MATERIAUX
!     -------------------------------------------------
    nomres(1) = 'ALPHA'
    nomres(2) = 'MU'
    nomres(3) = 'PORO'
    nomres(4) = 'KAPA'
    nomres(5) = 'LAMBDA'
    nomres(6) = 'M'
    nomres(7) = 'PRES_CRIT'
    nomres(8) = 'KCAM'
    nomres(9) = 'PTRAC'
!
    nompar(1) = 'TEMP'
    valpam(1) = tm
!
! - Get elastic properties
!
    call rcvala(imate, ' ', 'ELAS', 1, nompar, &
                valpam, 1, nomres(1), valres(1), icodre(1), &
                0)
!
! - Compute thermic dilation
!
    if (.not. isnan(tp) .and. isnan(tm)) then
        if (isnan(tref) .or. icodre(1) .ne. 0) then
            call utmess('F', 'COMPOR5_42')
        else
            coef = valres(1)*(tp-tref)-valres(1)*(tm-tref)
        end if
    else
        valres(1) = 0.d0
        coef = 0.d0
    end if
!
! - Get material properties
!
    call rcvala(imate, ' ', 'CAM_CLAY ', 1, nompar, &
                valpam, 8, nomres(2), valres(2), icodre(2), &
                2)
    mu = valres(2)
    poro = valres(3)
    kapa = valres(4)
    lambda = valres(5)
    m = valres(6)
    prescr = valres(7)
    kcam = valres(8)
    ptrac = valres(9)
    deuxmu = deux*mu
    e0 = poro/(1.d0-poro)
    xk0 = (1.d0+e0)/kapa
    xk = (1.d0+e0)/(lambda-kapa)
    if (kcam .ne. zero .and. kcam .le. (-xk0*ptrac)) then
        call utmess('F', 'COMPOR1_42')
    end if
!
!     -- 3 CALCUL DE DEPSMO ET DEPSDV :
!     --------------------------------
    if (cplan) then
        call utmess('F', 'ALGORITH6_63')
    end if
    depsmo = 0.d0
    do k = 1, ndimsi
        depsth(k) = deps(k)
    end do
    do k = 1, 3
        depsth(k) = depsth(k)-coef
        depsmo = depsmo+depsth(k)
    end do
    depsmo = -depsmo
    do k = 1, ndimsi
        depsdv(k) = depsth(k)+depsmo/3.d0*kron(k)
    end do
!
!     -- 4 CALCUL DE SIGMMO, SIGMDV, SIGEL,SIMOEL,SIELEQ, SIEQM :
!     -------------------------------------------------------------
    sigmmo = 0.d0
    do k = 1, 3
        sigmmo = sigmmo+sigm(k)
    end do
    sigmmo = -sigmmo/3.d0
    if (sigmmo .lt. ptrac) then
        call utmess('F', 'ALGORITH6_64')
    end if
    sieleq = 0.d0
    sieqm = 0.d0
    do k = 1, ndimsi
        sigmdv(k) = sigm(k)+sigmmo*kron(k)
        sieqm = sieqm+sigmdv(k)**2
        sigel(k) = sigmdv(k)+deuxmu*depsdv(k)
        sieleq = sieleq+sigel(k)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
    sieqm = sqrt(1.5d0*sieqm)
!
    if (((xk0*depsmo) .gt. epxmax)) then
        retcom = 1
        goto 999
    end if
    simoel = sigmmo*exp(xk0*depsmo)+kcam/xk0*(exp(xk0*depsmo)-un)
! ---- INITIALISATION A T=0
    if (pcrm(1) .eq. 0.d0) then
!
        pcrm(1) = prescr
        pcrm(3) = simoel
        pcrm(4) = sieleq
        pcrm(5) = 0.d0
        pcrm(6) = 0.d0
        pcrm(7) = e0
!
! ---- ON VERIFIE LA COHERENCE DES DONNEES MECA DE DEPART
        nu = (trois*((un+e0)*sigmmo+kapa*kcam)-deuxmu*kapa)/ &
             (six*((un+e0)*sigmmo+kapa*kcam)+deuxmu*kapa)
        young = deuxmu*(un+nu)
!
        if ((young .le. zero) .or. (nu .le. zero) .or. (nu .gt. unsde)) then
            call tecael(iadzi, iazk24)
            nomail = zk24(iazk24-1+3) (1:8)
            call utmess('F', 'COMPOR1_3', sk=nomail)
        end if
!
    end if
!
!
!     -- 5 CALCUL DU CRITERE :
!     ----------------------
    fonc = sieleq**2+m*m*(simoel-ptrac)**2-2.d0*m*m*(simoel-ptrac)*pcrm(1)
!
!     -- 6  TEST DE PLASTIFICATION ET CALCUL DE PCRP SIGP, SIGPDV :
!     ------------------------------------------------------------
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        if (fonc .le. 0.d0) then
!      -- TRAITEMENT DE L'ELASTICITE
            pcrpLoc(1) = pcrm(1)
            pcrpLoc(2) = 0.d0
            do k = 1, ndimsi
                sigpdv(k) = sigel(k)
                sigp(k) = sigel(k)-simoel*kron(k)
            end do
!
            pcrpLoc(3) = simoel
            pcrpLoc(4) = sieleq
            pcrpLoc(5) = 0.d0
            pcrpLoc(6) = 0.d0
!
            pcrpLoc(7) = pcrm(7)
            if (pcrm(3) .ne. zero) then
                if (pcrpLoc(3)/pcrm(3) .gt. zero) then
                    pcrpLoc(7) = pcrm(7)-kapa*log(pcrpLoc(3)/pcrm(3))
                end if
            end if
!
        else
!     -- PLASTIFICATION : CALCUL DE LA DEFORMATION
!     -- VOLUMIQUE PLASTIQUE : DEPPMO
            pcrpLoc(2) = 1.d0
            seuil = m**2*(pcrm(1)-ptrac)**2
!
            xinf = 0.d0
!     -- RECHERCHE DE LA BORNE SUP
!
            if (abs((simoel-pcrm(1)-ptrac)/seuil) .lt. ptit) then
                v0 = 0.d0
                goto 100
            end if
!
            if (abs(xk0*simoel+kcam+xk*pcrm(1)) .lt. ptit) then
                if ((-deux*(simoel-pcrm(1)-ptrac)/(xk0*simoel+kcam-xk*ptrac)) .lt. zero) then
                    xsup = 1.d0/(xk+xk0)*log(abs(simoel-ptrac)/pcrm(1))
                else
!       RESULTAT D UN DEVELOPPEMENT LIMITE D ORDRE 2
                    if ((simoel-ptrac) .gt. pcrm(1)) then
                        xsup = +sqrt((-deux*(simoel-pcrm(1)-ptrac)/(xk0*simoel+kcam-xk*ptrac)))
                    else
                        xsup = -sqrt((-deux*(simoel-pcrm(1)-ptrac)/(xk0*simoel+kcam-xk*ptrac)))
                    end if
                end if
            else
!       RESULTAT D UN DEVELOPPEMENT LIMITE D ORDRE 1
                xsup = (simoel-pcrm(1)-ptrac)/(xk0*simoel+kcam+xk*pcrm(1))
            end if
!
!     --RESOLUTION AVEC LA METHODE DE NEWTON ENTRE LES BORNES
            v0 = xinf
            if (((-xk0*v0) .gt. epxmax) .or. ((xk*v0) .gt. epxmax)) then
                retcom = 1
                goto 999
            end if
!
            f1 = (simoel+kcam/xk0)*exp(-xk0*v0)-kcam/xk0-ptrac
            f2 = (simoel+kcam/xk0)*exp(-xk0*v0)-kcam/xk0-ptrac-2.d0*pcrm(1)*exp(xk*v0)
            f3 = (simoel+kcam/xk0)*exp(-xk0*v0)-kcam/xk0-ptrac-pcrm(1)*exp(xk*v0)
            f4 = (1.d0+3.d0*deuxmu*v0/2.d0/m/m/f3)
!
            f = sieleq**2+m**2*f4**2*f1*f2
!
!
            f1p = -(xk0*simoel+kcam)*exp(-xk0*v0)
            f2p = -(xk0*simoel+kcam)*exp(-xk0*v0)-2.d0*xk*pcrm(1)*exp(xk*v0)
            f3p = -(xk0*simoel+kcam)*exp(-xk0*v0)-xk*pcrm(1)*exp(xk*v0)
            f4p = 3.d0*deuxmu/2.d0/(m**2)*(f3-v0*f3p)/f3/f3
!
            fp = m**2*f4**2*(f1p*f2+f1*f2p)+2.d0*m**2*f4*f4p*f1*f2
!
!
            do iter = 1, nint(carcri(1))
!
!     --CRITERE DE CONVERGENCE
                if ((abs(f)/seuil) .le. carcri(3)) goto 100
!
!     --CONSTRUCTION DU NOUVEL ITERE
                v0 = v0-f/fp
                if (xsup .gt. 0.d0) then
                    if (v0 .le. xinf .or. v0 .ge. xsup) v0 = (xinf+xsup)/2
                else
                    if (v0 .le. xsup .or. v0 .ge. xinf) v0 = (xinf+xsup)/2
                end if
!
!     --CALCUL DE LA FONCTION EN V0 ET DE SA DERIVEE
!
                if (((-xk0*v0) .gt. epxmax) .or. ((xk*v0) .gt. epxmax)) then
                    retcom = 1
                    goto 999
                end if
                f1 = (simoel+kcam/xk0)*exp(-xk0*v0)-kcam/xk0-ptrac
                f2 = (simoel+kcam/xk0)*exp(-xk0*v0)-kcam/xk0-ptrac-2.d0*pcrm(1)*exp(xk*v0)
                f3 = (simoel+kcam/xk0)*exp(-xk0*v0)-kcam/xk0-ptrac-pcrm(1)*exp(xk*v0)
                f4 = (1.d0+3.d0*deuxmu*v0/2.d0/m/m/f3)
!
                f = sieleq**2+m**2*f4**2*f1*f2
!
                if (f .gt. zero) signf = 1
                if (f .lt. zero) signf = -1
!
                f1p = -(xk0*simoel+kcam)*exp(-xk0*v0)
                f2p = -(xk0*simoel+kcam)*exp(-xk0*v0)-2.d0*xk*pcrm(1)*exp(xk*v0)
                f3p = -(xk0*simoel+kcam)*exp(-xk0*v0)-xk*pcrm(1)*exp(xk*v0)
                f4p = 3.d0*deuxmu/2.d0/(m**2)*(f3-v0*f3p)/f3/f3
!
                fp = m**2*f4**2*(f1p*f2+f1*f2p)+2.d0*m**2*f4*f4p*f1*f2
!
!
                if (((-xk0*xinf) .gt. epxmax) .or. ((xk*xinf) .gt. epxmax)) then
                    retcom = 1
                    goto 999
                end if
!
                fxi1 = (simoel+kcam/xk0)*exp(-xk0*xinf)-kcam/xk0-ptrac
                fxi2 = (simoel+kcam/xk0)*exp(-xk0*xinf)-kcam/xk0-ptrac-2.d0*pcrm(1)*exp(xk*xinf)
                fxi3 = (simoel+kcam/xk0)*exp(-xk0*xinf)-kcam/xk0-ptrac-pcrm(1)*exp(xk*xinf)
                fxi4 = (1.d0+3.d0*deuxmu*xinf/2.d0/m/m/fxi3)
!
                fxi = sieleq**2+m**2*fxi4**2*fxi1*fxi2
!
                if (fxi .gt. zero) signfi = 1
                if (fxi .lt. zero) signfi = -1
!
                if ((signf*signfi) .lt. zero) xsup = v0
                if ((signf*signfi) .gt. zero) xinf = v0
!
            end do
            retcom = 1
            goto 999
100         continue
            deppmo = v0
!
!     -- REACTUALISATION DE LA VARIABLE INTERNE
            if (((xk*deppmo) .gt. epxmax) .or. (xk0*(depsmo-deppmo) .gt. epxmax)) then
                retcom = 1
                goto 999
            end if
!
!
            pcrpLoc(1) = pcrm(1)*exp(xk*deppmo)
!     -- REACTUALISATION DES CONTRAINTES
            sigpmo = (sigmmo+kcam/xk0)*exp(xk0*(depsmo-deppmo))-kcam/xk0
            call r8inir(6, 0.d0, sigpdv, 1)
            do k = 1, ndimsi
                sigpdv(k) = sigel(k)/ &
                            (1.d0+(3.d0*deuxmu/2.d0*deppmo)/(m*m*(sigpmo-pcrpLoc(1)-ptrac)))
                sigp(k) = sigpdv(k)-sigpmo*kron(k)
            end do
!
!
! ---- V(3) CONTRAINTE VOLUMIQUE
            pcrpLoc(3) = sigpmo
!
! ---- V(4) CONTRAINTE EQUIVALENTE
            sieqp = 0.0d0
            do k = 1, ndimsi
                sieqp = sieqp+sigpdv(k)**2.d0
            end do
            pcrpLoc(4) = sqrt(1.5d0*sieqp)
!
! ---- V(5) DEFORMATION PLASTIQUE VOLUMIQUE
            pcrpLoc(5) = pcrm(5)+deppmo
!
! ---- V(6) DEFORMATION PLASTIQUE EQUIVALENTE
            depseq = 0.0d0
            do k = 1, ndimsi
                depseq = depseq+depsdv(k)*depsdv(k)
            end do
            depseq = sqrt(2.d0/3.d0*depseq)
            pcrpLoc(6) = pcrm(6)+depseq
!
! ---- V(7) :: INDICE DES VIDES
            pcrpLoc(7) = pcrm(7)
            if (pcrm(3) .ne. zero .and. pcrm(1) .ne. zero) then
                if (pcrpLoc(3)/pcrm(3) .gt. zero .and. pcrpLoc(1)/pcrm(1) .gt. zero) then
                    pcrpLoc(7) = pcrm(7)- &
                                 kapa*log(pcrpLoc(3)/pcrm(3))- &
                                 (lambda-kapa)*log(pcrpLoc(1)/pcrm(1))
                end if
            end if
        end if
!
    end if
!
!
!     -- 7 CALCUL DE L'OPERATEUR TANGENT :
!     --------------------------------
    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
!
        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
            if (pcrm(2) .eq. 0.d0) then
                matr = 0
            else
                matr = 1
            end if
        end if
        if (option(1:9) .eq. 'FULL_MECA') then
            if (pcrpLoc(2) .eq. 1.d0) then
                matr = 2
            else
                matr = 0
            end if
        end if
!      INITIALISATION DE L'OPERATEUR TANGENT
!     ---------------------------------------
        dsidep(:, :) = 0.d0
!
!     -- 7.1 CALCUL DE DSIDEP(6,6)-ELASTIQUE:
!     ---------------------------------------
        if (matr .eq. 0) then
            do k = 1, 3
                do l = 1, 3
                    dsidep(k, l) = xk0*simoel+kcam-deuxmu/3.d0
                end do
            end do
            do k = 1, ndimsi
                dsidep(k, k) = dsidep(k, k)+deuxmu
            end do
        end if
!
!     -- 7.2 CALCUL DE DSIDEP(6,6)-EN VITESSE :
!     ---------------------------------------
        if (matr .eq. 1) then
!
            call r8inir(6*6, 0.d0, dsidep, 1)
!     -- 7.2.1 CALCUL DU MODULE ELASTOPLASTIQUE H
!
            valm = 0.d0
            do i = 1, ndimsi
                valm = valm+sigmdv(i)**2
            end do
!
            h = 4.d0*m**4*(sigmmo-ptrac)*(sigmmo-ptrac-pcrm(1))* &
                (xk0*(sigmmo-ptrac-pcrm(1))+xk*pcrm(1))+deuxmu*9.d0*valm
!
!
!     -- 7.2.2 CALCUL D'UN TERME INTERMEDIAIRE
            do k = 1, 3
                a(k) = 0.d0
            end do
            do k = 1, 3
                a(k) = -deux*xk0*m*m*(sigmmo-ptrac)*(sigmmo-ptrac-pcrm(1))*kron(k)+ &
                       3.d0*deuxmu*sigmdv(k)
            end do
            call r8inir(3, 0.d0, aa, 1)
            do k = 4, ndimsi
                aa(k) = 3.d0*deuxmu*sigmdv(k)
            end do
!
!     -- 7.2.3 CALCUL DES TERMES DE DSIDEP
            call r8inir(ndimsi*ndimsi, 0.d0, dsidep, 1)
            do k = 1, 3
                do l = 1, 3
                    dsidep(k, l) = xk0*(sigmmo-ptrac)-deuxmu/3.d0-a(k)*a(l)/h
                end do
            end do
            do k = 1, 3
                do l = 4, ndimsi
                    dsidep(k, l) = -a(k)*aa(l)
                    dsidep(k, l) = dsidep(k, l)/h
                    dsidep(l, k) = dsidep(k, l)
                end do
            end do
            do k = 4, ndimsi
                do l = 4, ndimsi
                    dsidep(k, l) = -aa(k)*aa(l)
                    dsidep(k, l) = dsidep(k, l)/h
                end do
            end do
            do k = 1, ndimsi
                dsidep(k, k) = deuxmu+dsidep(k, k)
            end do
!
        end if
!
        if (matr .eq. 2) then
            call r8inir(6*6, 0.d0, dsidep, 1)
!
!     -- 7.2.1 CALCUL DU MODULE ELASTOPLASTIQUE H
!
            valp = 0.d0
            do i = 1, ndimsi
                valp = valp+sigpdv(i)**2
            end do
!
            h = 4.d0*m**4*(sigpmo-ptrac)* &
                (sigpmo-ptrac-pcrpLoc(1))* &
                (xk0*(sigpmo-ptrac-pcrpLoc(1))+xk*pcrpLoc(1))+deuxmu*9.d0*valp
!
!     -- 7.2.2 CALCUL D'UN TERME INTERMEDIAIRE
            call r8inir(3, 0.d0, a, 1)
            call r8inir(3, 0.d0, aa, 1)
!
            do k = 1, 3
                a(k) = -deux*xk0*m*m*(sigpmo-ptrac)*(sigpmo-ptrac-pcrpLoc(1))*kron(k)+ &
                       3.d0*deuxmu*sigpdv(k)
            end do
!
            do k = 4, ndimsi
                aa(k) = 3.d0*deuxmu*sigpdv(k)
            end do
!
!     -- 7.2.3 CALCUL DES TERMES DE DSIDEP
            call r8inir(ndimsi*ndimsi, 0.d0, dsidep, 1)
            do k = 1, 3
                do l = 1, 3
                    dsidep(k, l) = xk0*(sigpmo-ptrac)-deuxmu/3.d0-a(k)*a(l)/h
                end do
            end do
            do k = 1, 3
                do l = 4, ndimsi
                    dsidep(k, l) = -a(k)*aa(l)
                    dsidep(k, l) = dsidep(k, l)/h
                    dsidep(l, k) = dsidep(k, l)
                end do
            end do
            do k = 4, ndimsi
                do l = 4, ndimsi
                    dsidep(k, l) = -aa(k)*aa(l)
                    dsidep(k, l) = dsidep(k, l)/h
                end do
            end do
            do k = 1, ndimsi
                dsidep(k, k) = deuxmu+dsidep(k, k)
            end do
!
        end if
!     -- 7.3 CALCUL DE DSIDEP(6,6)-MATRICE COHERENTE :
!     ----------------------------------------------
        if (matr .eq. 3) then
            sieqp = 0.0d0
            do k = 1, ndimsi
                sieqp = sieqp+sigpdv(k)**2
            end do
            sieqp = sqrt(1.5d0*sieqp)
            diff2 = abs((pcrpLoc(1)-sigpmo)/pcrpLoc(1))
            if (diff2 .lt. carcri(3)) then
!
!     -- 7.3.1 OPERATEUR TANGENT COHERENT AU POINT CRITIQUE
!     -- TRAITEMENT DE LA PARTIE DEVIATORIQUE
!     -- CALCUL DE Q+
!     -- CALCUL DU TENSEUR HH QUI MULTIPLIE LA DEFORMATION
                call r8inir(6*6, 0.d0, ses, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        ses(k, l) = 1.d0/2.d0*(sigpdv(k)*sigel(l)+sigel(k)*sigpdv(l))
                    end do
                end do
                call r8inir(6*6, 0.d0, hh, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        hh(k, l) = -deuxmu*3.d0*ses(k, l)/2.d0/sieleq/sieqp
                    end do
                end do
                do k = 1, ndimsi
                    hh(k, k) = deuxmu+hh(k, k)
                end do
                if (ndim .eq. 2) then
                    hh(5, 5) = 1.d0
                    hh(6, 6) = 1.d0
                end if
!     -- INVERSE DE HH
                call r8inir(6*6, 0.d0, hhm, 1)
                do k = 1, 6
                    hhm(k, k) = 1.d0
                end do
                call mgauss('NFWP', hh, hhm, 6, 6, 6, rbid, iret)
!
!     -- CALCUL DU TENSEUR GG QUI MULTIPLIE LA CONTRAINTE
                call r8inir(6*6, 0.d0, gg, 1)
                call r8inir(6*6, 0.d0, sps, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        sps(k, l) = sigpdv(k)*sigpdv(l)
                    end do
                end do
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        gg(k, l) = -3.d0*sieleq*sps(k, l)/2.d0/sieqp**3
                    end do
                end do
                do k = 1, ndimsi
                    gg(k, k) = sieleq/sieqp+gg(k, k)
                end do
!     --  MATRICE DE PROJECTION SUR L'ESPACE DES CONTRAINTES
!     -- DEVIATORIQUES
                call r8inir(6*6, 0.d0, v, 1)
                do k = 1, 3
                    do l = 1, 3
                        v(k, l) = -1.d0/3.d0
                        v(l, k) = v(k, l)
                    end do
                end do
                do k = 1, ndimsi
                    v(k, k) = v(k, k)+1.d0
                end do
!     --  PRODUIT DE LA MATRICE DE PROJECTION SUR L'ESPACE
!     --  DES CONTRAINTES DEVIATORIQUES PAR GG
                call r8inir(6*6, 0.d0, d1g, 1)
                call promat(v, 6, ndimsi, ndimsi, gg, &
                            6, ndimsi, ndimsi, d1g)
!     -- PRODUIT DU RESULTAT PAR L'INVERSE DE HH
                call r8inir(6*6, 0.d0, d1ghhm, 1)
                call promat(d1g, 6, ndimsi, ndimsi, hhm, &
                            6, ndimsi, ndimsi, d1ghhm)
!
!     -- 7.3.2 TRAITEMENT DE LA PARTIE HYDROSTATIQUE
!     --  PRODUIT DE LA MATRICE DE PROJECTION SUR L'ESPACE
!     --  DES CONTRAINTES DEVIATORIQUES PAR LA MATRICE IDENTITE
!     --  D'ORDRE 2
                call r8inir(6*6, 0.d0, id2, 1)
                do k = 1, 3
                    do l = 1, 3
                        id2(k, l) = -1.d0/3.d0/xk0/sigpmo
                    end do
                end do
!     -- SOMME DES TERMES DEVIATORIQUE ET HYDROSTATIQUE
                call r8inir(6*6, 0.d0, devhyd, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        devhyd(k, l) = d1ghhm(k, l)/deuxmu+id2(k, l)
                    end do
                end do
                if (ndim .eq. 2) then
                    devhyd(5, 5) = 1.d0
                    devhyd(6, 6) = 1.d0
                end if
!     -- INVERSE DE LA SOMME DES TERMES DEVIATORIQUE ET HYDROSTATIQUE
                call r8inir(6*6, 0.d0, devhym, 1)
                do k = 1, 6
                    devhym(k, k) = 1.d0
                end do
                call mgauss('NFWP', devhyd, devhym, 6, 6, 6, rbid, iret)
!     -- TERMES DE L'OPERATEUR TANGENT
                call r8inir(6*6, 0.d0, dsidep, 1)
                do k = 1, 6
                    do l = 1, 6
                        dsidep(k, l) = devhym(k, l)
                    end do
                end do
            else
!
!      ---7.4 OPERATEUR TANGENT COHERENT CAS GENERAL
!      -- CALCUL DES INCREMENTS DE P ET DE S
                deltap = sigpmo-sigmmo
                call r8inir(6, 0.d0, deltas, 1)
                do k = 1, ndimsi
                    deltas(k) = sigpdv(k)-sigmdv(k)
                end do
!
!     --  CALCUL DE VECTEURS INTERMEDIAIRES
                spards = 0.d0
                do k = 1, ndimsi
                    spards = spards+deltas(k)*sigpdv(k)
                end do
                call r8inir(6, 0.d0, tplus, 1)
                do k = 1, ndimsi
                    tplus(k) = sigpdv(k)+deltas(k)
                end do
!
!      -- 7.4.1 TERMES NECESSAIRES A LA PARTIE DEVIATORIQUE
                hp = 4.d0*m**4*xk*sigpmo*pcrpLoc(1)*(sigpmo-pcrpLoc(1))
!
                xc = 9.d0*spards/hp
                xd = 6.d0*m*m*(sigpmo-pcrpLoc(1))*deltap/hp
                xv = 3.d0*spards+2.d0*m**2*(sigpmo-pcrpLoc(1))*deltap
                xlam = xv/hp
                xa = (4.d0*xlam*xk*m**4*sigpmo*(sigpmo-2.d0*pcrpLoc(1))+ &
                      2.d0*m**2*deltap)*m**2*(sigpmo-pcrpLoc(1))/ &
                     (m**2*xlam+(1.d0/2.d0/xk/pcrpLoc(1)))
                xi = 2.d0*m**2*(sigpmo-pcrpLoc(1))- &
                     2.d0*m**4*xlam*(sigpmo-pcrpLoc(1))/((1.d0/2.d0/xk/pcrpLoc(1))+m**2*xlam)
                rap = xi/(hp+xa)
!
!     -- CALCUL DE LA MATRICE CC-SYMETRISATION DE TPLUS.I
!
                call r8inir(6*6, 0.d0, cc, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        cc(k, l) = (tplus(k)*kron(l)+kron(k)*tplus(l))/2.d0
                    end do
                end do
!
!     -- CALCUL DES TERMES D'UNE MATRICE INTERMEDIAIRE C
!
                call r8inir(6*6, 0.d0, c, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        c(k, l) = 9.d0/2.d0/(hp+xa)*(sigpdv(k)*tplus(l)+tplus(k)*sigpdv(l))
                    end do
                end do
                do k = 1, ndimsi
                    c(k, k) = c(k, k)+1.d0/deuxmu+xc+xd
                end do
!
!     -- ASSEMBLAGE DES TERMES POUR LA PARTIE DEVIATORIQUE
                call r8inir(6*6, 0.d0, ee, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        ee(k, l) = c(k, l)-rap*cc(k, l)
                    end do
                end do
!
!      -- TERMES NECESSAIRES A LA PARTIE HYDROSTATIQUE
                xu = 2.d0*m**2*xk*pcrpLoc(1)
                xg = xlam*xu/(1.d0+xlam*xu)
                xh = xu*(sigpmo-pcrpLoc(1))/(1.d0+xlam*xu)
                xe = 1.d0+xh*2.d0*m**2*deltap/hp+xh*4.d0*xk*m**4*sigpmo* &
                     (sigpmo-2.d0*pcrpLoc(1))*xv/hp/hp
                xf = (2.d0*m**2*(sigpmo-pcrpLoc(1))+2.d0*m**2*deltap-xg*2.d0*m**2*deltap)/hp- &
                     4.d0*xk*m**4*xv/hp/hp*((2.d0*sigpmo-pcrpLoc(1))*pcrpLoc(1)+ &
                                            xg*sigpmo*(sigpmo-2.d0*pcrpLoc(1)))
                ct = (1.d0+ &
                      2.d0*m**2*xk0*sigpmo*( &
                      xlam-xg*xlam-xlam*xf*xh/xe+xf/xe*(sigpmo-pcrpLoc(1)) &
                      ))/(xk0*sigpmo)
!
!     --  VECTEUR INTERMEDIAIRE
                call r8inir(6, 0.d0, fv, 1)
                do k = 1, ndimsi
                    fv(k) = 3.d0*xf/xe*sigpdv(k)-ct*kron(k)/3.d0
                end do
!     -- SYMMETRISATION DEFV ET SA PROJECTION SUR L'ESPACE
!     -- DES CONTRAINTES HYDROSTATIQUES
                call r8inir(6*6, 0.d0, ffi, 1)
                do k = 1, 3
                    do l = 1, 3
                        ffi(k, l) = -(1.d0/3.d0)*(fv(k)+fv(l))/2.d0
                    end do
                end do
                do k = 1, 3
                    do l = 4, ndimsi
                        ffi(k, l) = -(1.d0/3.d0)*fv(l)/2.d0
                        ffi(l, k) = ffi(k, l)
                    end do
                end do
!     --  MATRICE DE PROJECTION SUR L'ESPACE DES CONTRAINTES
!     -- DEVIATORIQUES
                call r8inir(6*6, 0.d0, v, 1)
                do k = 1, 3
                    do l = 1, 3
                        v(k, l) = -1.d0/3.d0
                        v(l, k) = v(k, l)
                    end do
                end do
                do k = 1, ndimsi
                    v(k, k) = v(k, k)+1.d0
                end do
!     -- PROJECTION DE EE SUR L'ESPACE DES CONTRAINTES
!     -- DEVIATORIQUES
                call r8inir(6*6, 0.d0, s, 1)
                call promat(ee, 6, ndimsi, ndimsi, v, &
                            6, ndimsi, ndimsi, s)
!
!     -- COMBINAISON DES DEUX PARTIES DEVIATORIQUE ET
!     -- HYDROSTATIQUE
                call r8inir(6*6, 0.d0, t, 1)
                do k = 1, ndimsi
                    do l = 1, ndimsi
                        t(k, l) = s(k, l)+ffi(k, l)
                    end do
                end do
                if (ndim .eq. 2) then
                    t(5, 5) = 1.d0
                    t(6, 6) = 1.d0
                end if
!     -- INVERSE DE LA MATRICE T
                call r8inir(6*6, 0.d0, vv, 1)
                do k = 1, 6
                    vv(k, k) = 1.d0
                end do
                call mgauss('NFWP', t, vv, 6, 6, &
                            6, rbid, iret)
!     --  7.3.3 CALCUL DES TERMES DSIDEP L'OPERATEUR TANGENT
                call r8inir(6*6, 0.d0, dsidep, 1)
                do k = 1, 6
                    do l = 1, 6
                        dsidep(k, l) = vv(k, l)
                    end do
                end do
!
            end if
        end if
    end if
! ======================================================================
999 continue

    if (lVari) then
        pcrp = pcrpLoc
    end if
! =====================================================================
end subroutine
