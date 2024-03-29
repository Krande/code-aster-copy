! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine huresi(mod, nmat, mater, indi, deps, &
                  nr, yd, yf, nvi, vind, &
                  r, iret)

! aslint: disable=W1306
    implicit none
!     ----------------------------------------------------------------
!     CALCUL DU VECTEUR RESIDU DU SYSTEME NL
!     ----------------------------------------------------------------
!     IN   MOD    :  TYPE DE MODELISATION
!          NMAT   :  DIMENSION TABLEAU DONNEES MATERIAU
!          MATER  :  DONNEES MATERIAU
!          INDI   :  MECANISMES POTENTIELLEMENT ACTIFS
!          DEPS   :  INCREMENT DEFORMATION
!          YD     :  VECTEUR SOLUTION A T
!          YF     :  VECTEUR SOLUTION A T+DT?
!          VIND   :  VARIABLES INTERNES A T
!          NVI    :  NOMBRE DE VARIABLES INTERNES
!          NR     :  DIMENSION MAXIMALE DE YD
!     OUT  R      :  VECTEUR RESIDU DU SYSTEME NL A RESOUDRE
!          IRET   :  CODE RETOUR (>0 -> PB)
!     ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/hujddd.h"
#include "asterfort/hujksi.h"
#include "asterfort/hujpic.h"
#include "asterfort/hujprc.h"
#include "asterfort/hujprj.h"
#include "asterfort/hujpxd.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
    integer :: indi(7), nr, nvi, iret, nmat
    character(len=8) :: mod
    real(kind=8) :: mater(nmat, 2), deps(6), yd(nr), yf(nr), vind(nvi), r(*)
    integer :: ndt, ndi, i, nbmeca, nbmect, k, j, kk
    real(kind=8) :: le(6), ccond, levp, lr(4), pref, lf(7), hooknl(6, 6)
    real(kind=8) :: depse(6), cde(6), coef0, i1f, n, hook(6, 6)
    real(kind=8) :: beta, d, b, phi, angdil, pco, acyc, amon, ccyc, cmon
    real(kind=8) :: m, mdil, coef, ptrac, zero, degr, trois
    real(kind=8) :: e, nu, al, demu, la, un, e1, e2, e3, nu12, nu23, nu13, deux
    real(kind=8) :: g1, g2, g3, nu21, nu31, nu32, denom, depsp(6), dlambd(7)
    real(kind=8) :: psi(42), sigf(6), sigdc(9), ad(7), ksi(7), q(7), p(7)
    real(kind=8) :: rc(7), dpsids(6, 6), sigd(3), th(2), prod, tole1
    real(kind=8) :: yft(nr), ydt(nr), mul, ps, pk, rtrac, dpsi
    real(kind=8) :: epsvp, pc, matert(22, 2)
    aster_logical :: prox(4), proxc(4)
!
    parameter(ndt=6)
    parameter(ndi=3)
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    parameter(degr=0.0174532925199d0)
    parameter(tole1=1.d-6)
!     ----------------------------------------------------------------
! --- REDIMENSIONNEMENT DE YD ET YF POUR S'ADAPTER A HUJJID
! --- COPIE A PARTIR DU TRAITEMENT DE HUJMID
    ydt(1:nr) = yd(1:nr)
    yft(1:nr) = yf(1:nr)
!
    do i = 1, 22
        matert(i, 1) = mater(i, 1)
        matert(i, 2) = mater(i, 2)
    end do
!
    do i = 1, 6
        ydt(i) = yd(i)*mater(1, 1)
        yft(i) = yf(i)*mater(1, 1)
    end do
!
    nbmeca = 0
    nbmect = 0
    do k = 1, 7
        if (indi(k) .gt. 0) then
            nbmect = nbmect+1
            if (indi(k) .le. 8) nbmeca = nbmeca+1
        end if
        dlambd(k) = zero
        ad(k) = zero
        ksi(k) = zero
        q(k) = zero
        p(k) = zero
    end do
!
    do i = 1, nbmeca
        ydt(ndt+1+i) = yd(ndt+1+i)*mater(1, 1)/abs(mater(8, 2))
        yft(ndt+1+i) = yf(ndt+1+i)*mater(1, 1)/abs(mater(8, 2))
    end do
! ====================================================================
! --- PROPRIETES HUJEUX MATERIAU -------------------------------------
! ====================================================================
    n = mater(1, 2)
    beta = mater(2, 2)
    d = mater(3, 2)
    b = mater(4, 2)
    phi = mater(5, 2)
    angdil = mater(6, 2)
    pco = mater(7, 2)
    pref = mater(8, 2)
    acyc = mater(9, 2)
    amon = mater(10, 2)
    ccyc = deux*mater(11, 2)
    cmon = mater(12, 2)
    m = sin(degr*phi)
    mdil = sin(degr*angdil)
    coef = mater(20, 2)
    ptrac = mater(21, 2)
    rtrac = abs(pref*1.d-6)
!
! ====================================================================
! --- OPERATEURS ELASTICITE LINEAIRES---------------------------------
! ====================================================================
    hook(:, :) = zero
!
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN') then
!
        if (mater(17, 1) .eq. un) then
!
            e = mater(1, 1)
            nu = mater(2, 1)
            al = e*(un-nu)/(un+nu)/(un-deux*nu)
            demu = e/(un+nu)
            la = e*nu/(un+nu)/(un-deux*nu)
!
            do i = 1, ndi
                do j = 1, ndi
                    if (i .eq. j) hook(i, j) = al
                    if (i .ne. j) hook(i, j) = la
                end do
            end do
            do i = ndi+1, ndt
                hook(i, i) = demu
            end do
!
        else if (mater(17, 1) .eq. deux) then
!
            e1 = mater(1, 1)
            e2 = mater(2, 1)
            e3 = mater(3, 1)
            nu12 = mater(4, 1)
            nu13 = mater(5, 1)
            nu23 = mater(6, 1)
            g1 = mater(7, 1)
            g2 = mater(8, 1)
            g3 = mater(9, 1)
            nu21 = mater(13, 1)
            nu31 = mater(14, 1)
            nu32 = mater(15, 1)
            denom = mater(16, 1)
!
            hook(1, 1) = (un-nu23*nu32)*e1/denom
            hook(1, 2) = (nu21+nu31*nu23)*e1/denom
            hook(1, 3) = (nu31+nu21*nu32)*e1/denom
            hook(2, 2) = (un-nu13*nu31)*e2/denom
            hook(2, 3) = (nu32+nu31*nu12)*e2/denom
            hook(3, 3) = (un-nu21*nu12)*e3/denom
            hook(2, 1) = hook(1, 2)
            hook(3, 1) = hook(1, 3)
            hook(3, 2) = hook(2, 3)
            hook(4, 4) = g1*2.d0
            hook(5, 5) = g2*2.d0
            hook(6, 6) = g3*2.d0
!
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
        call utmess('F', 'COMPOR1_4')
    end if
! ====================================================================
! --- OPERATEUR ELASTICITE NON LINEAIRE ------------------------------
! ====================================================================
    i1f = trace(ndi, yft)/trois
    if ((i1f/pref) .lt. tole1) i1f = tole1*pref
!
    coef0 = (i1f/pref)**n
    do i = 1, ndt
        do j = 1, ndt
            hooknl(i, j) = coef0*hook(i, j)
        end do
    end do
! ====================================================================
! --- AUTRES GRANDEURS UTILES ----------------------------------------
! ====================================================================
    do i = 1, 4
        prox(i) = .false.
        proxc(i) = .false.
    end do
!
    do i = 1, ndt
        sigf(i) = yft(i)
        psi(i) = zero
        psi(ndt+i) = zero
        psi(2*ndt+i) = zero
        psi(3*ndt+i) = zero
        psi(4*ndt+i) = zero
        psi(5*ndt+i) = zero
        psi(6*ndt+i) = zero
    end do
!
    do i = 1, 9
        sigdc(i) = zero
    end do
!
    do k = 1, nbmect
        kk = indi(k)
!
!        IF(YFT(NDT+1+NBMECA+K).GT.ZERO)THEN
        dlambd(k) = yft(ndt+1+nbmeca+k)
!        ENDIF
!
        if (kk .le. 8) rc(k) = yft(ndt+1+k)
!
        call hujddd('PSI   ', indi(k), matert, indi, yft, &
                    vind, psi((k-1)*ndt+1), dpsids, iret)
        if (iret .eq. 1) goto 998
!
        if (indi(k) .lt. 4) then
!
            call hujprj(indi(k), sigf, sigd, p(k), q(k))
            if (p(k) .ge. ptrac) goto 997
            call hujksi('KSI   ', matert, rc(k), ksi(k), iret)
            if (iret .eq. 1) goto 998
            ad(k) = acyc+ksi(k)*(amon-acyc)
!
        else if (indi(k) .eq. 4) then
!
            ksi(k) = un
            p(k) = trace(ndi, yft)/trois
!
        else if ((indi(k) .lt. 8) .and. (indi(k) .gt. 4)) then
!
            call hujprc(k, indi(k)-4, sigf, vind, matert, &
                        yft, p(k), q(k), sigdc(3*k-2))
            if (p(k) .ge. ptrac) goto 997
            call hujksi('KSI   ', matert, rc(k), ksi(k), iret)
            if (iret .eq. 1) goto 998
            ad(k) = deux*(acyc+ksi(k)*(amon-acyc))
!
            th(1) = vind(4*indi(k)-9)
            th(2) = vind(4*indi(k)-8)
            prod = sigdc(3*k-2)*th(1)+sigdc(3*k)*th(2)/deux
!
            if ((-q(k)/pref .lt. tole1) .or. ((un+prod/q(k)) .lt. tole1)) then
                kk = kk-4
                call hujpxd(indi(k), matert, sigf, vind, prox(kk), &
                            proxc(kk))
            else
                ad(k) = (acyc+ksi(k)*(amon-acyc))*(un+prod/q(k))
            end if
!
        else if (indi(k) .eq. 8) then
!
            ksi(k) = un
            call hujpic(k, indi(k), sigf, vind, matert, &
                        yft, p(k))
!
        else if ((indi(k) .gt. 8) .and. (indi(k) .lt. 12)) then
            cycle
!
        else
            call utmess('F', 'COMPOR1_8')
        end if
!
    end do
!
    epsvp = yft(ndt+1)
    pc = pco*exp(-beta*epsvp)
    cmon = cmon*pc/pref
    ccyc = ccyc*pc/pref
!
! --- CONDITIONNEMENT DE LA MATRICE JACOBIENNE
    ccond = mater(1, 1)
! ====================================================================
! ---- CALCUL DE CDE = C*DEPSE ---------------------------------------
! ====================================================================
    do i = 1, ndt
        depsp(i) = zero
    end do
!
    do k = 1, nbmect
        kk = (k-1)*ndt
        do i = 1, ndt
            depsp(i) = depsp(i)+dlambd(k)*psi(kk+i)
        end do
    end do
!
    do i = 1, ndt
        depse(i) = deps(i)-depsp(i)
    end do
!
    cde(1:ndt) = matmul(hooknl(1:ndt, 1:ndt), depse(1:ndt))
! ====================================================================
! --- CALCUL DE LE (6) -----------------------------------------------
! ====================================================================
    do i = 1, ndt
        le(i) = yft(i)-ydt(i)-cde(i)
    end do
!
!
! ====================================================================
! --- CALCUL DE LEVP (1X1) -------------------------------------------
! ====================================================================
    levp = yft(ndt+1)-ydt(ndt+1)
!
    if (nbmeca .eq. 0) goto 190
!
    do k = 1, nbmect
!
        kk = indi(k)
        pk = p(k)-ptrac
!
        if (kk .lt. 4) then
!
            if ((p(k)/pref) .gt. tole1) then
                dpsi = mdil+q(k)/p(k)
            else
                dpsi = mdil+q(k)
            end if
            levp = levp+coef*dlambd(k)*ksi(k)*dpsi
!
        else if (kk .eq. 4) then
!
            levp = levp+dlambd(k)
!
        else if ((kk .gt. 4) .and. (kk .lt. 8)) then
!
            call hujprj(kk-4, sigf, sigd, coef0, mul)
            ps = 2.d0*sigd(1)*sigdc(3*k-2)+sigd(3)*sigdc(3)
!
            if ((p(k)/pref) .gt. tole1) then
                if ((-q(k)/pref) .gt. tole1) then
                    dpsi = mdil+ps/2.d0/p(k)/q(k)
                else
                    dpsi = mdil
                end if
            else
                if ((-q(k)/pref) .gt. tole1) then
                    dpsi = mdil+ps/2.d-6/pref/q(k)
                else
                    dpsi = mdil
                end if
            end if
!
            levp = levp+coef*dlambd(k)*ksi(k)*dpsi
!
        else if (kk .eq. 8) then
!
            if (vind(22) .gt. zero) then
                levp = levp-dlambd(k)
            else
                levp = levp+dlambd(k)
            end if
!
        end if
    end do
!
! ====================================================================
! --- CALCUL DE LR (NBMECX1) -----------------------------------------
! ====================================================================
    do k = 1, 4
        lr(k) = zero
    end do
!
    if (nbmeca .eq. 0) goto 190
    do k = 1, nbmeca
        kk = indi(k)
        if (kk .lt. 4) then
            lr(k) = yft(ndt+1+k)-ydt(ndt+1+k)-dlambd(k)/ad(k)*(un-rc(k))**deux
        else if (kk .eq. 4) then
            lr(k) = yft(ndt+1+k)-ydt(ndt+1+k)-dlambd(k)/cmon*(un-rc(k))**deux
!
        else if ((kk .gt. 4) .and. (kk .lt. 8)) then
            th(1) = vind(4*indi(k)-9)
            th(2) = vind(4*indi(k)-8)
            prod = sigdc(3*k-2)*th(1)+sigdc(3*k)*th(2)/deux
!
            if ((-q(k)/pref .lt. tole1) .or. ((un+prod/q(k)) .lt. tole1)) then
                ad(k) = (acyc+ksi(k)*(amon-acyc))
            else
                ad(k) = (acyc+ksi(k)*(amon-acyc))*(un+prod/q(k))
            end if
            lr(k) = yft(ndt+1+k)-ydt(ndt+1+k)-dlambd(k)/ad(k)*(un-rc(k))**deux
        else if (kk .eq. 8) then
            lr(k) = yft(ndt+1+k)-ydt(ndt+1+k)-dlambd(k)/ccyc*(un-rc(k))**deux
!
        end if
    end do
!
190 continue
!
! ====================================================================
! --- CALCUL DE LF (NBMECX1) -----------------------------------------
! ====================================================================
    do k = 1, 7
        lf(k) = zero
    end do
!
    do k = 1, nbmect
        kk = indi(k)
        pk = p(k)-ptrac
        if (kk .lt. 4) then
            lf(k) = q(k)+m*pk*rc(k)*(un-b*log(pk/pc))
        else if (kk .eq. 4) then
            lf(k) = abs(p(k))+rc(k)*d*pc
        else if ((kk .gt. 4) .and. (kk .lt. 8)) then
            lf(k) = q(k)+m*pk*rc(k)*(un-b*log(pk/pc))
        else if (kk .eq. 8) then
            lf(k) = abs(p(k))+rc(k)*d*pc
        else if (kk .gt. 8) then
            call hujprj(kk-8, yft, sigd, pk, ps)
            lf(k) = pk+deux*rtrac-ptrac
        end if
    end do
!
! ====================================================================
! --- ASSEMBLAGE DE R ------------------------------------------------
! ====================================================================
    do i = 1, ndt
        r(i) = -le(i)/ccond
    end do
!
    r(ndt+1) = -levp
!
    if (nbmeca .eq. 0) goto 240
!
    do k = 1, nbmeca
        r(ndt+1+k) = -lr(k)/ccond*abs(pref)
        r(ndt+1+nbmeca+k) = -lf(k)/ccond
    end do
!
240 continue
!
    if (nbmeca .lt. nbmect) then
        do k = 1, nbmect
            if (indi(k) .gt. 8) then
                r(ndt+1+nbmeca+k) = -lf(k)/ccond
            end if
        end do
    end if
!
    do i = ndt+nbmeca+nbmect+2, 18
        r(i) = zero
    end do
!
    goto 999
! ====================================================================
! --- TRAITEMENT DES CAS OU IRET DIFFERENT DE ZERO -------------------
! ====================================================================
997 continue
    iret = 3
!
998 continue
!
    do i = 1, 3
        if (prox(i)) then
            vind(i+4) = mater(18, 2)
            vind(23+i) = un
            vind(27+i) = zero
            vind(4*i+5) = zero
            vind(4*i+6) = zero
            vind(4*i+7) = zero
            vind(4*i+8) = zero
            vind(5*i+31) = zero
            vind(5*i+32) = zero
            vind(5*i+33) = zero
            vind(5*i+34) = zero
            vind(5*i+35) = mater(18, 2)
            iret = 2
        else if (proxc(i)) then
            vind(27+i) = zero
            iret = 2
        end if
    end do
!
999 continue
!
end subroutine
