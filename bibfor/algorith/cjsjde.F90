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

subroutine cjsjde(mod, mater, epsd, deps, yd, &
                  yf, gd, r, signe, drdy)
! aslint: disable=W1501
    implicit none
!     INTEGRATION PLASTIQUE (MECANISME DEVIATOIRE SEUL) DE LA LOI CJS
!
!     RESOLUTION PAR METHODE DE NEWTON       DRDY(DYI) DDYI = - R(DYI)
!
!     CALCUL DU SECOND MEMBRE : - R(DYI)
!     CALCUL DU JACOBIEN      : DRDY(DYI)
!                        Y =  (SIG     ,  R      ,  X      ,  LAMBD  )
!                        R = -( LE     ,  LCR    ,  LCX    ,  FD     )
!                     DRDY =  ( DLEDS  ,  DLEDR  ,  DLEDX  ,  DLEDL  )
!                             ( DLRDS  ,  DLRDR  ,  DLRDX  ,  DLRDL  )
!                             ( DLXDS  ,  DLXDR  ,  DLXDX  ,  DLXDL  )
!                             ( DFDDS  ,  DFDDR  ,  DFDDX  ,  DFDDL  )
!     ------------------------------------------------------------------
!     IN   MOD      :  MODELISATION
!          MATER    :  COEFFICIENTS MATERIAU A T+DT
!          EPSD     :  DEFORMATION A T
!          DEPS     :  INCREMENT DE DEFORMATION
!          YD       :  VARIABLES A T = (SIGD, VIND, LAMBDD)
!          YF       :  VARIABLES A T+DT = (SIGF, VINF, LAMBDF)
!     VAR  GD       :  TENSEUR DE LA LOI D ECOULEMENT PLASTIQUE DEV.
!     OUT  R        :  SECOND MEMBRE
!          SIGNE    :  SIGNE DE S:DEPSDP
!          DRDY     :  JACOBIEN
! ======================================================================
#include "asterfort/calcq.h"
#include "asterfort/cjsdtd.h"
#include "asterfort/cjsqco.h"
#include "asterfort/cjst.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcicma.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nmod, i, j, k, codret
    parameter(nmod=14)
    real(kind=8) :: epsd(6), deps(6), depse(6), depsdp(6)
    real(kind=8) :: dsignl(6), dsigl(6), sigf(6)
    real(kind=8) :: yd(nmod), yf(nmod), r(nmod), drdy(nmod, nmod)
    real(kind=8) :: mater(14, 2), n, rm, rc, a, b, c, pco, pa, pc
    real(kind=8) :: beta, betapr, gamma, mucjs
    real(kind=8) :: le(6), lcr, lcx(6), fd
    real(kind=8) :: hooknl(6, 6), hook(6, 6), kooh(6, 6)
!        REAL*8       KOOHNL(6,6)
    real(kind=8) :: i1f, e, nu, al, la, mu, unpnue, unsure, mnuse
    real(kind=8) :: rf, gr, dgrdr, dgrds(6), dhdq(6)
    real(kind=8) :: xf(6), xii, gx(6), dgxdx(6, 6), dgxds(6, 6)
    real(kind=8) :: epsv, phi, phio, rr, cosa, cosdif
    real(kind=8) :: truc, signe, dlambd
    real(kind=8) :: s(6), sii, siic, hts, dets, cos3ts
    real(kind=8) :: siirel, qiirel
    real(kind=8) :: q(6), qii, htq, detq, cos3tq
    real(kind=8) :: tangs, tangq, tetas, tetaq
    real(kind=8) :: qq(6), qqii, norm(6), vectan(6), gd(6), gdd(6)
    real(kind=8) :: coef0, coef1, coef4, coef5, coef6
    real(kind=8) :: coef7, coef8, coef9, coef10
    real(kind=8) :: coef13, coef14, coef15, coef16, coef17, coef18
    real(kind=8) :: coef19, coef20, coef21, coef22
    real(kind=8) :: dgdds(6, 6), dgddr(6), dgddx(6, 6)
    real(kind=8) :: dleds(6, 6), dledr(6), dledx(6, 6), dledl(6)
    real(kind=8) :: dlrds(6), dlrdr, dlrdx(6), dlrdl
    real(kind=8) :: dlxds(6, 6), dlxdr(6), dlxdx(6, 6), dlxdl(6)
    real(kind=8) :: dfdds(6), dfddr, dfddx(6), dfddl
    real(kind=8) :: dvds(6, 6)
    real(kind=8) :: dsds(6, 6), dssds(6), ds2ds(6), ds2cds(6)
    real(kind=8) :: dqqds(6, 6), dqiids(6), dqq2ds(6), dqqdx(6, 6)
    real(kind=8) :: dcads(6), dcfds(6), drrds(6), dphods(6), dphids(6)
    real(kind=8) :: dhtsds(6), dhtqds(6)
    real(kind=8) :: t(6), td(6), dtddq(6, 6), ts(6), dhds(6)
    real(kind=8) :: dqqdq(6, 6), dqds(6, 6)
    real(kind=8) :: d2fdsq(6, 6), d2fds2(6, 6), d2fdsx(6, 6)
    real(kind=8) :: terme1, terme2, terme3
    real(kind=8) :: prod0, prod1, prod2, prod3, prod4, prod5, prod6
    real(kind=8) :: mun5, zero, un, d12, deux, trois, cinq, six
    real(kind=8) :: kron(6), iden6(6, 6), quatre
    real(kind=8) :: epssig, pref, qinit
    character(len=8) :: mod
! ======================================================================
    parameter(mun5=-1.5d0)
    parameter(d12=0.5d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    parameter(quatre=4.d0)
    parameter(cinq=5.d0)
    parameter(six=6.d0)
    parameter(epssig=1.d-8)
! ======================================================================
    common/tdim/ndt, ndi
    data iden6/un, zero, zero, zero, zero, zero,&
         &                   zero, un, zero, zero, zero, zero,&
         &                   zero, zero, un, zero, zero, zero,&
         &                   zero, zero, zero, un, zero, zero,&
         &                   zero, zero, zero, zero, un, zero,&
         &                   zero, zero, zero, zero, zero, un/
    data kron/un, un, un, zero, zero, zero/
! ======================================================================
    call jemarq()
! ======================================================================
! --- ATTENTION : NE PAS CONFONDRE LA VARIABLE INTERNE R ---------------
! --- ---------   AVEC LE SYSTEME NON LINEAIRE GLOBAL NOTE AUSSI R -----
! ======================================================================
! --- PROPRIETES CJS MATERIAU ------------------------------------------
! ======================================================================
    beta = mater(1, 2)
    rm = mater(2, 2)
    n = mater(3, 2)
    rc = mater(5, 2)
    a = mater(6, 2)
    b = mater(7, 2)
    c = mater(8, 2)
    gamma = mater(9, 2)
    mucjs = mater(10, 2)
    pco = mater(11, 2)
    pa = mater(12, 2)
    qinit = mater(13, 2)
! ======================================================================
! --- PREMIER INVARIANT ET AUTRES GRANDEURS UTILES ---------------------
! ======================================================================
    i1f = trace(ndi, yf)
    if ((i1f+qinit) .eq. 0.d0) then
        i1f = -qinit+1.d-12*pa
        pref = abs(pa)
    else
        pref = abs(i1f+qinit)
    end if
!
    rf = yf(ndt+1)
!
    do i = 1, ndt
        xf(i) = yf(ndt+1+i)
        gdd(i) = gd(i)
    end do
!
    dlambd = yf(2*ndt+2)-yd(2*ndt+2)
!
    do i = 1, ndt
        sigf(i) = yf(i)
    end do
! ======================================================================
! --- OPERATEURS DE RIGIDITE ET DE SOUPLESSE (LINEAIRES OU NON LINEA.) -
! ======================================================================
! --- OPERATEURS LINEAIRES ---------------------------------------------
! ======================================================================
    hook(:, :) = zero
    kooh(:, :) = zero
    e = mater(1, 1)
    nu = mater(2, 1)
    al = e*(un-nu)/(un+nu)/(un-deux*nu)
    la = nu*e/(un+nu)/(un-deux*nu)
    mu = e*d12/(un+nu)
! ======================================================================
    unpnue = (un+nu)/e
    unsure = un/e
    mnuse = -nu/e
! ======================================================================
! --- 3D/DP/AX ---------------------------------------------------------
! ======================================================================
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        do i = 1, ndi
            do j = 1, ndi
                if (i .eq. j) then
                    hook(i, j) = al
                    kooh(i, j) = unsure
                end if
                if (i .ne. j) then
                    hook(i, j) = la
                    kooh(i, j) = mnuse
                end if
            end do
        end do
        do i = ndi+1, ndt
            do j = ndi+1, ndt
                if (i .eq. j) then
                    hook(i, j) = deux*mu
                    kooh(i, j) = unpnue
                end if
            end do
        end do
! ======================================================================
! --- CP/1D ------------------------------------------------------------
! ======================================================================
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'ALGORITH2_15')
    end if
! ======================================================================
! --- OPERATEURS NON LINEAIRE ------------------------------------------
! ======================================================================
    coef0 = ((i1f+qinit)/trois/pa)**n
!
    do i = 1, ndt
        do j = 1, ndt
            hooknl(i, j) = coef0*hook(i, j)
        end do
    end do
! ======================================================================
! --- LOIS D ECROUISSAGE : GR ET GX ------------------------------------
! ======================================================================
! --- ECROUISSAGE ISOTROPE ---------------------------------------------
! ======================================================================
    coef1 = ((i1f+qinit)/trois/pa)**mun5
    gr = -a*(un-rf/rm)**deux*(i1f+qinit)*coef1
! ======================================================================
! --- ECROUISSAGE CINEMATIQUE ------------------------------------------
! ======================================================================
    call cjsqco(gamma, sigf, xf, pref, epssig, &
                i1f, s, sii, siirel, cos3ts, &
                hts, dets, q, qii, qiirel, &
                cos3tq, htq, detq)
!
    call calcq(q, gamma, pref, epssig, qq, &
               codret)
    qqii = norm2(qq(1:ndt))
!
    xii = norm2(xf(1:ndt))
!
    epsv = zero
    do i = 1, ndi
        epsv = epsv+epsd(i)+deps(i)
    end do
!
    pc = pco*exp(-c*epsv)
!
    if (xii .le. epssig) then
        phi = un
    else if (siirel .le. epssig) then
        cosa = un
        cosdif = un
        rr = rc+mucjs*max(zero, log(trois*pc/(i1f+qinit)))
        phio = cosa/(rr-hts/htq*rm*cosdif)
        phi = phio*hts*qqii
    else
        cosa = (qii*qii-sii*sii-i1f*i1f*xii*xii)/(deux*sii*i1f*xii)
!
        tangs = sqrt(un-cos3ts*cos3ts)/cos3ts
        tangq = sqrt(un-cos3tq*cos3tq)/cos3tq
        tetas = atan2(tangs, 1.d0)/trois
        tetaq = atan2(tangq, 1.d0)/trois
        cosdif = cos(tetas-tetaq)
!
        rr = rc+mucjs*max(zero, log(trois*pc/(i1f+qinit)))
        phio = cosa/(rr-hts/htq*rm*cosdif)
        phi = phio*hts*qqii
    end if
!
    do i = 1, ndt
        gx(i) = (i1f+qinit)/b*(qq(i)+phi*xf(i))*coef1
    end do
! ======================================================================
! --- LOI D ECOULEMENT DU MECANISME PLASTIQUE DEVIATOIRE : GD ----------
! ======================================================================
! --- LOI D ECOULEMENT -------------------------------------------------
! ======================================================================
!
    truc = dot_product(qq(1:ndt), xf(1:ndt))-rf
!
    do i = 1, ndt
        dfdds(i) = qq(i)-truc*kron(i)
    end do
! ======================================================================
! --- CALCUL DE L INCREMENT DE DEFORMATION PLASTIQUE DEV. EN -----------
! --- UTILISANT LE TENSEUR GD DE L ITERATION DE NEWTON PRECEDENTE ------
! --- ET LA VALEUR DE DLAMBD -------------------------------------------
! ======================================================================
    do i = 1, ndt
        depsdp(i) = dlambd*gdd(i)
    end do
!
    truc = dot_product(s(1:ndt), depsdp(1:ndt))
    if (truc .ge. zero) then
        signe = un
    else
        signe = -un
    end if
!
    siic = -rc*(i1f+qinit)/hts
    betapr = beta*(sii/siic-un)*signe
    coef4 = betapr/sii
    coef5 = un/sqrt(betapr*betapr+trois)
!
    do i = 1, ndt
        vectan(i) = coef4*s(i)+kron(i)
        norm(i) = coef5*vectan(i)
    end do
!
    prod0 = dot_product(dfdds(1:ndt), norm(1:ndt))
    do i = 1, ndt
        gd(i) = dfdds(i)-prod0*norm(i)
    end do
! ======================================================================
! --- CALCULS PRELIMIAIRES DE DERIVEES ---------------------------------
! ======================================================================
! --- DERIVEE DE S PAR RAPPORT A SIG : DSDS ----------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dsds(i, j) = iden6(i, j)-kron(i)*kron(j)/trois
        end do
    end do
! ======================================================================
! --- DERIVEE DE Q PAR RAPPORT A SIG : DQDS ----------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dqds(i, j) = iden6(i, j)-kron(j)*(kron(i)/trois+xf(i))
        end do
    end do
! ======================================================================
! --- EXPRESSION DE TS, T, TD ET DE SA DERIVEE PAR RAPPORT A Q, DTDDQ --
! ======================================================================
    call cjst(s, ts)
    call cjst(q, t)
    call lcdevi(t, td)
    call cjsdtd(mod, q, dtddq)
! ======================================================================
! --- DERIVEE DE HTQ PAR RAPPORT A Q : DHDQ ----------------------------
! --- DERIVEE DE HTS PAR RAPPORT A S : DHDS ----------------------------
! ======================================================================
    coef6 = sqrt(trois/deux)*gamma/htq**cinq/qii**trois
    coef7 = -gamma*cos3tq/(deux*htq**cinq*qii**deux)
    coef8 = sqrt(trois/deux)*gamma/hts**cinq/sii**trois
    coef9 = -gamma*cos3ts/(deux*hts**cinq*sii**deux)
    do i = 1, ndt
        dhdq(i) = coef6*t(i)+coef7*q(i)
        dhds(i) = coef8*ts(i)+coef9*s(i)
    end do
! ======================================================================
! --- DERIVEE DE QQ PAR RAPPORT A Q : DQQDQ ----------------------------
! ======================================================================
    coef6 = un+gamma/deux*cos3tq
    coef7 = sqrt(54.0d0)*gamma/(six*qii*qii)
    coef8 = sqrt(54.0d0)*gamma/deux/qii**quatre/htq**cinq
    coef9 = sqrt(54.0d0)*gamma/six/qii**deux/htq**cinq
!
    do i = 1, ndt
        do j = 1, ndt
            dqqdq(i, j) = -cinq/htq**six*(coef6*q(i)/qii+coef7*td(i))* &
                          dhdq(j)+coef6/htq**cinq*(iden6(i, j)/qii-q(i)*q(j)/qii** &
                                            trois)+coef8*q(i)*(t(j)-3*q(j)*detq/qii**deux)+coef9*( &
                          dtddq(i, j)-deux*td(i)*q(j)/qii**deux)
        end do
    end do
! ======================================================================
! --- DERIVEE DE QQ PAR RAPPORT A SIG : DQQDS --------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dqqds(i, j) = zero
            do k = 1, ndt
                dqqds(i, j) = dqqds(i, j)+dqqdq(i, k)*dqds(k, j)
            end do
        end do
    end do
! ======================================================================
! --- DERIVEE DE DFDDS PAR RAPPORT A Q : D2FDSQ ------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            prod1 = zero
            do k = 1, ndt
                prod1 = prod1+dqqdq(k, j)*xf(k)
            end do
            d2fdsq(i, j) = dqqdq(i, j)-prod1*kron(i)
        end do
    end do
! ======================================================================
! --- DERIVEE DE DFDDS PAR RAPPORT A SIG : D2FDS2 ----------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            d2fds2(i, j) = zero
            do k = 1, ndt
                d2fds2(i, j) = d2fds2(i, j)+d2fdsq(i, k)*dqds(k, j)
            end do
        end do
    end do
! ======================================================================
! --- DERIVEE DE SII PAR RAPPORT A SIG : DS2DS -------------------------
! ======================================================================
    do i = 1, ndt
        ds2ds(i) = zero
        do j = 1, ndt
            ds2ds(i) = ds2ds(i)+s(j)*dsds(i, j)
        end do
        ds2ds(i) = ds2ds(i)/sii
    end do
! ======================================================================
! --- DERIVEE DE SIIC PAR RAPPORT A SIG : DS2CDS -----------------------
! ======================================================================
    coef9 = -rc/hts
    coef10 = rc*(i1f+qinit)/hts**deux
    do i = 1, ndt
        prod2 = zero
        do j = 1, ndt
            prod2 = prod2+dhds(j)*dsds(j, i)
        end do
        ds2cds(i) = coef9*kron(i)+coef10*prod2
    end do
! ======================================================================
! --- DERIVEE DU RAPPORT (SII/SIIC) PAR RAPPORT A SIG : DSSDS ----------
! ======================================================================
    do i = 1, ndt
        dssds(i) = ds2ds(i)/siic-ds2cds(i)*sii/siic/siic
    end do
! ======================================================================
! --- DERIVEE DU VECTEUR TANGENT A LA SURFACE POTENTIELLE, VECTAN, -----
! --- PAR RAPPORT A SIG : DVDS -----------------------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dvds(i, j) = betapr/sii*dsds(i, j)+signe*beta*s(i)*(ds2ds(j)/sii/sii-ds2cds(j)&
                 &/siic/siic)
        end do
    end do
! ======================================================================
! --- DERIVEE DE QII PAR RAPPORT A SIG : DQIIDS ------------------------
! ======================================================================
    do i = 1, ndt
        dqiids(i) = zero
        do j = 1, ndt
            dqiids(i) = dqiids(i)+q(j)/qii*dqds(j, i)
        end do
    end do
! ======================================================================
! --- DERIVEE DE QQII PAR RAPPORT A SIG : DQQ2DS -----------------------
! ======================================================================
    do i = 1, ndt
        dqq2ds(i) = zero
        do j = 1, ndt
            dqq2ds(i) = dqq2ds(i)+qq(j)/qqii*dqqds(j, i)
        end do
    end do
! ======================================================================
! --- DERIVEE DE HTQ PAR RAPPORT A SIG : DHTQDS ------------------------
! --- DERIVEE DE HTS PAR RAPPORT A SIG : DHTSDS ------------------------
! ======================================================================
    do i = 1, ndt
        dhtqds(i) = zero
        dhtsds(i) = zero
        do j = 1, ndt
            dhtqds(i) = dhtqds(i)+dhdq(j)*dqds(j, i)
            dhtsds(i) = dhtsds(i)+dhds(j)*dsds(j, i)
        end do
    end do
! ======================================================================
! --- DERIVEE DE PHI PAR RAPPORT A SIG : DPHIDS ------------------------
! --- EN FONCTION DES DIFFERENTES VALEURS DE SII ET XII ----------------
! ======================================================================
! --- INITIALISATION : -------------------------------------------------
! ======================================================================
    do i = 1, ndt
        dphids(i) = zero
    end do
! ======================================================================
! --- 1ER CAS : XII = 0  ON NE FAIT RIEN -------------------------------
! ======================================================================
    if (xii .gt. epssig) then
! ======================================================================
! --- 2EME CAS : SII = 0 -----------------------------------------------
! ======================================================================
        if (siirel .le. epssig) then
! ======================================================================
! --- DERIVEE DE RR PAR RAPPORT A SIG : DRRDS --------------------------
! ======================================================================
            do i = 1, ndt
                drrds(i) = -mucjs/(i1f+qinit)*kron(i)
            end do
! ======================================================================
! --- DERIVEE DE PHIO PAR RAPPORT A SIG : DPHODS -----------------------
! ======================================================================
            do i = 1, ndt
                dphods(i) = -cosa/(rr-hts/htq*rm*cosdif)**deux*(drrds(i)-rm*cosdif/htq*&
                     &dhtsds(i)+hts/htq/htq*rm*cosdif*dhtqds(i))
            end do
! ======================================================================
! --- DERIVEE DE PHI PAR RAPPORT A SIG : DPHIDS ------------------------
! ======================================================================
            do i = 1, ndt
                dphids(i) = hts*qqii*dphods(i)+phio*qqii*dhtsds(i)+phio*hts*dqq2ds(i)
            end do
! ======================================================================
! --- 3EME CAS : SII ET XII NON NULS -----------------------------------
! ======================================================================
        else
! ======================================================================
! --- DERIVEE DE COSA PAR RAPPORT A SIG : DCADS ------------------------
! ======================================================================
            do i = 1, ndt
                dcads(i) = ( &
                     qii*dqiids(i)-i1f*xii*xii*kron(i)-sii*ds2ds(i))/sii/i1f/xii-d12*( &
                     &qii*qii-sii*sii-i1f*i1f*xii*xii)/(sii*i1f*xii)**deux*(i1f*xii*d&
                     &s2ds(i)+sii*xii*kron(i) &
                     )
            end do
! ======================================================================
! --- DERIVEE DE COSDIF PAR RAPPORT A SIG : DCFDS ----------------------
! ======================================================================
            coef17 = sqrt(54.0d0)/sii**trois
            coef18 = sqrt(54.0d0)/qii**trois
            coef19 = trois*dets/sii**deux
            coef20 = trois*detq/qii**deux
            coef21 = sqrt(un-cos3ts*cos3ts)
            coef22 = sqrt(un-cos3tq*cos3tq)
!
            do i = 1, ndt
                dcfds(i) = sin(tetas-tetaq)/trois*(coef21*(coef17*(ts(i)-coef19*s(i)&
                     & ))-coef22*(coef18*(t(i)-coef20*q(i))))
            end do
! ======================================================================
! --- DERIVEE DE RR PAR RAPPORT A SIG : DRRDS --------------------------
! ======================================================================
            do i = 1, ndt
                drrds(i) = -mucjs/(i1f+qinit)*kron(i)
            end do
! ======================================================================
! --- DERIVEE DE PHIO PAR RAPPORT A SIG : DPHODS -----------------------
! ======================================================================
            do i = 1, ndt
                dphods(i) = dcads(i)/(rr-hts/htq*rm*cosdif)-cosa/(rr-hts/htq*rm*cosdif&
                     & )**deux*(drrds(i)-rm*cosdif/htq*dhtsds(i)+hts/htq/htq*rm*cos&
                     &dif*dhtqds(i)-hts/htq*rm*dcfds(i))
            end do
! ======================================================================
! --- DERIVEE DE PHI PAR RAPPORT A SIG : DPHIDS ------------------------
! ======================================================================
            do i = 1, ndt
                dphids(i) = hts*qqii*dphods(i)+phio*qqii*dhtsds(i)+phio*hts*dqq2ds(i)
            end do
        end if
    end if
! ======================================================================
! --- DERIVEE DE QQ PAR RAPPORT A X : DQQDX ----------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dqqdx(i, j) = zero
            do k = 1, ndt
                dqqdx(i, j) = dqqdx(i, j)-i1f*dqqdq(i, k)*iden6(k, j)
            end do
        end do
    end do
! ======================================================================
! --- DERIVEE DE DFDDS PAR RAPPORT A X : D2FDSX ------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            prod3 = zero
            do k = 1, ndt
                prod3 = prod3+dqqdx(k, j)*xf(k)+qq(k)*iden6(k, j)
            end do
            d2fdsx(i, j) = dqqdx(i, j)-prod3*kron(i)
        end do
    end do
! ======================================================================
! --- DERIVEES DES LOIS D ECROUISSAGE GR ET GX -------------------------
! --- PAR RAPPORT A R, X ET SIG :  DGRDS, DGRDR, DGXDS, DGXDX ----------
! --- ET DERIVEES DE LA LOI D ECOULEMENT GD : DGDDS, DGDDR, DGDDX ------
! ======================================================================
! --- DERIVEES DE LA LOI D ECROUISSAGE ISOTROPE ------------------------
! ======================================================================
    dgrdr = deux*a/rm*(un-rf/rm)*(i1f+qinit)*coef1
    do i = 1, ndt
        dgrds(i) = a/deux*(un-rf/rm)**deux*coef1*kron(i)
    end do
! ======================================================================
! --- DERIVEES DE LA LOI D ECROUISSAGE CINEMATIQUE ---------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dgxds(i, j) = -coef1/deux/b*( &
                 qq(i)+phi*xf(i))*kron(j)+(i1f+qinit)*coef1/b*(dqqds(i, j)+xf(&
                 &i)*dphids(j) &
                 )
            dgxdx(i, j) = (i1f+qinit)*coef1/b*(dqqdx(i, j)+phi*iden6(i, j))
        end do
    end do
! ======================================================================
! --- DERIVEES DE LA LOI D ECOULEMENT ----------------------------------
! ======================================================================
    coef13 = coef5*coef5*deux*beta*beta*(sii/siic-un)
    do i = 1, ndt
        do j = 1, ndt
            prod4 = zero
            prod5 = zero
            prod6 = zero
            do k = 1, ndt
                prod4 = prod4+d2fds2(k, j)*norm(k)
                prod5 = prod5+dfdds(k)*dvds(k, j)
                prod6 = prod6+d2fdsx(k, j)*norm(k)
            end do
!
            dgdds(i, j) = d2fds2(i, j)-prod4*norm(i)-coef5*prod5*norm(i)-coef5*prod0 &
                         &*dvds(i, j)+coef13*prod0*norm(i)*dssds(j)
            dgddx(i, j) = d2fdsx(i, j)-prod6*norm(i)
        end do
    end do
!
    coef14 = betapr*betapr/(betapr*betapr+trois)
    coef15 = -trois*betapr/sii/(betapr*betapr+trois)
!
    do i = 1, ndt
        dgddr(i) = coef14*kron(i)+coef15*s(i)
    end do
! ======================================================================
! --- LOI D ETAT : LE --------------------------------------------------
! --- ET DERIVEE DE LA LOI D ETAT : DLEDS, DLEDR, DLEDX, DLEDL ---------
! ======================================================================
! --- LOI D ETAT -------------------------------------------------------
! ======================================================================
    do i = 1, ndt
        depse(i) = deps(i)-dlambd*gd(i)
    end do
!
    dsignl(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), depse(1:ndt))
!
    do i = 1, ndt
        le(i) = yf(i)-yd(i)-dsignl(i)
    end do
! ======================================================================
! --- DERIVEE DE LA LOI D ETAT -----------------------------------------
! ======================================================================
    coef16 = n/trois/pa*(trois*pa/(i1f+qinit))**(un-n)
    dsigl(1:ndt) = matmul(hook(1:ndt, 1:ndt), depse(1:ndt))
!
    dleds(:, :) = zero
!
    do i = 1, ndt
        do j = 1, ndt
            terme1 = zero
            terme2 = zero
            do k = 1, ndt
                terme1 = terme1+hooknl(i, k)*dgdds(k, j)
                terme2 = terme2+hooknl(i, k)*dgddx(k, j)
            end do
            dleds(i, j) = iden6(i, j)-coef16*dsigl(i)*kron(j)+dlambd*terme1
            dledx(i, j) = dlambd*terme2
        end do
!
        terme3 = zero
        do k = 1, ndt
            terme3 = terme3+hooknl(i, k)*dgddr(k)
        end do
        dledr(i) = dlambd*terme3
    end do
    dledl(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), gd(1:ndt))
! ======================================================================
! - LOI D ECROUISSAGE DE R : LCR ---------------------------------------
! - ET DERIVEE DE LA LOI D ECROUISSAGE DE R: DLRDS, DLRDR, DLRDX, DLRDL-
! ======================================================================
! --- LOI D ECROUISSAGE DE R -------------------------------------------
! ======================================================================
    lcr = rf-yd(ndt+1)-dlambd*gr
! ======================================================================
! --- DERIVEE DE LA LOI D ECROUISSAGE DE R -----------------------------
! ======================================================================
    do i = 1, ndt
        dlrds(i) = -dlambd*dgrds(i)
        dlrdx(i) = zero
    end do
    dlrdr = un-dlambd*dgrdr
    dlrdl = -gr
! ======================================================================
! - LOI D ECROUISSAGE DE X : LCX ---------------------------------------
! - ET DERIVEE DE LA LOI D ECROUISSAGE DE X: DLXDS, DLXDR, DLXDX, DLXDL-
! ======================================================================
! --- LOI D ECROUISSAGE DE X -------------------------------------------
! ======================================================================
    do i = 1, ndt
        lcx(i) = xf(i)-yd(ndt+1+i)-dlambd*gx(i)
    end do
! ======================================================================
! --- DERIVEE DE LA LOI D ECROUISSAGE DE X -----------------------------
! ======================================================================
    do i = 1, ndt
        do j = 1, ndt
            dlxds(i, j) = -dlambd*dgxds(i, j)
            dlxdx(i, j) = iden6(i, j)-dlambd*dgxdx(i, j)
        end do
        dlxdr(i) = zero
        dlxdl(i) = -gx(i)
    end do
! ======================================================================
! --- SEUIL DEVIATOIRE : FD --------------------------------------------
! --- ET DERIVEES DU SEUIL DEVIATOIRE : DFDDS, DFDDR, DFDDX, DFDDL -----
! ======================================================================
! --- SEUIL DEVIATOIRE -------------------------------------------------
! ======================================================================
    fd = qii*htq+rf*(i1f+qinit)
! ======================================================================
! --- DERIVEES DU SEUIL DEVIATOIRE : -----------------------------------
! ======================================================================
! --- DFDDS : DEJA CALCULE CI-DESSUS -----------------------------------
! ======================================================================
    dfddr = (i1f+qinit)
!
    do i = 1, ndt
        dfddx(i) = zero
        do j = 1, ndt
            dfddx(i) = dfddx(i)-(i1f+qinit)*(qii*dhdq(j)+htq*q(j)/qii)*iden6(j, i)
        end do
    end do
!
    dfddl = zero
! ======================================================================
! --- ASSEMBLAGE DE R = -( LE     ,  LCR    ,  LCX    ,  FD     ) ------
! ======================================================================
! --- ASSEMBLAGE DE DRDY -----------------------------------------------
! ======================================================================
!
!                  DRDY =  ( DLEDS  ,  DLEDR  ,  DLEDX  ,  DLEDL  )
!                          ( DLRDS  ,  DLRDR  ,  DLRDX  ,  DLRDL  )
!                          ( DLXDS  ,  DLXDR  ,  DLXDX  ,  DLXDL  )
!                          ( DFDDS  ,  DFDDR  ,  DFDDX  ,  DFDDL  )
!
! ======================================================================
! --- ASSEMBLAGE DE R --------------------------------------------------
! ======================================================================
    do i = 1, ndt
        r(i) = -le(i)
        r(ndt+1+i) = -lcx(i)
    end do
    r(ndt+1) = -lcr
    r(2*ndt+2) = -fd
! ======================================================================
! --- ASSEMBLAGE DE DRDY -----------------------------------------------
! ======================================================================
    call lcicma(dleds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                1, 1)
    call lcicma(dledr, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+1)
    call lcicma(dledx, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                1, ndt+2)
    call lcicma(dledl, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                1, 2*ndt+2)
!
    call lcicma(dlrds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+1, 1)
    drdy(ndt+1, ndt+1) = dlrdr
    call lcicma(dlrdx, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+1, ndt+2)
    drdy(ndt+1, 2*ndt+2) = dlrdl
!
    call lcicma(dlxds, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+2, 1)
    call lcicma(dlxdr, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                ndt+2, ndt+1)
    call lcicma(dlxdx, 6, 6, ndt, ndt, &
                1, 1, drdy, nmod, nmod, &
                ndt+2, ndt+2)
    call lcicma(dlxdl, 6, 1, ndt, 1, &
                1, 1, drdy, nmod, nmod, &
                ndt+2, 2*ndt+2)
!
    call lcicma(dfdds, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                2*ndt+2, 1)
    drdy(2*ndt+2, ndt+1) = dfddr
    call lcicma(dfddx, 1, 6, 1, ndt, &
                1, 1, drdy, nmod, nmod, &
                2*ndt+2, ndt+2)
    drdy(2*ndt+2, 2*ndt+2) = dfddl
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
