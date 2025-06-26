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
subroutine nmtacr(mode, ndimsi, mat, sigel, vim, &
                  epm, dp, sp, xi, f, &
                  g, fds, gds, fdp, gdp, &
                  fdx, gdx, dpmax, sig, tang)
!
    implicit none
#include "blas/ddot.h"
    integer(kind=8) :: mode, ndimsi
    real(kind=8) :: mat(14), sigel(ndimsi), vim(9), epm(6), dp, sp, xi
    real(kind=8) :: f, g, fdp, fds, fdx, gdp, gds, gdx, dpmax, sig(ndimsi)
    real(kind=8) :: tang(6, 6)
!
!
! ----------------------------------------------------------------------
! TAHERI :  EVALUATION DES CRITERES ET DERIVEES
! ----------------------------------------------------------------------
! IN  MODE   SELECTION DES TERMES A CALCULER
!             0 -> DPMAX
!             1 -> F,G
!             2 -> F, G, FDS, GDS, FDP, GDP
!             3 -> F, G, FDS, GDS, FDX, GDX
!             4 -> SIG
!             5 -> TANG (PLASTICITE CLASSIQUE)
!             6 -> TANG (PLASTICITE A DEUX SURFACES)
! IN  NDIMSI DIMENSION DES TENSEURS
! IN  MAT    TABLEAU DES CONSTANTES MATERIAUX
! IN  SIGEL  DEVIATEUR DES CONTRAINTES ELASTIQUES
! IN  VIM    VARIABLES INTERNES EN T-
! IN  EPM    DEFORMATION PLASTIQUE EN T-
! IN  DP     INCREMENT DE DEFORMATION PLASTIQUE CUMULEE
! IN  SP     CONTRAINTE DE PIC (SIGMA P)
! IN  XI     PILOTAGE DE EPN
! OUT F      VALEUR DU CRITERE DE PLASTICITE
! OUT G      VALEUR DU CRITERE DE CONTRAINTE MAXIMALE
! OUT FDS    DERIVEE DE F PAR RAPPORT A SIGMA P
! OUT GDS    DERIVEE DE G PAR RAPPORT A SIGMA P
! OUT FDP    DERIVEE DE F PAR RAPPORT A P
! OUT GDP    DERIVEE DE G PAR RAPPORT A P
! OUT FDX    DERIVEE DE F PAR RAPPORT A XI
! OUT GDX    DERIVEE DE G PAR RAPPORT A XI
! OUT DPMAX  BORNE SUPERIEURE POUR L'INCREMENT DE DEFO. PLAST.
! OUT SIG    CORRECTION DE CONTRAINTE (DEUXMU * DEP)
! OUT TANG   VARIATION :  DEP/SIGEL * PROJDEV
! ----------------------------------------------------------------------
!
    integer(kind=8) :: n, i, j
    real(kind=8) :: epndx(6), epn(6)
    real(kind=8) :: c, d, a(6), b(6), se(6), s0(6), seeq
    real(kind=8) :: v(6), semax
    real(kind=8) :: p, seq, m(6), meq, l(6), leq, q
    real(kind=8) :: tmp, cdp, ddp, sedp(6), seeqdp, s0dp(6), seqdp
    real(kind=8) :: mdp(6), meqdp, qdp, ldp(6), leqdp
    real(kind=8) :: cds, dds, ads(6), seds(6), seeqds, s0ds(6), seqds
    real(kind=8) :: mds(6), meqds, qds, lds(6), leqds
    real(kind=8) :: adx(6), bdx(6), sedx(6), seeqdx, s0dx(6), seqdx
    real(kind=8) :: mdx(6), meqdx, ldx(6), leqdx, qdx
    real(kind=8) :: pbs0, seqde(6), meqde(6), qde(6), fde(6)
    real(kind=8) :: pas0, leqde(6), gde(6)
    real(kind=8) :: pde(6), det
    blas_int :: b_incx, b_incy, b_n
!
!
!
! -- CALCUL DES ELEMENTS COMMUNS
!
    p = vim(1)+dp
    tmp = exp(-mat(8)*p*(1.d0-sp/mat(11)))
    c = mat(10)+mat(9)*tmp
    d = 1.d0-mat(6)*tmp
    do n = 1, ndimsi
        epndx(n) = vim(2+n)-epm(n)
        b(n) = -xi*epndx(n)
        epn(n) = epm(n)-b(n)
        a(n) = (mat(11)-sp)*epm(n)+sp*b(n)
        se(n) = sigel(n)-c*a(n)
    end do
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    seeq = sqrt(1.5d0*ddot(b_n, se, b_incx, se, b_incy))
!
!
! -- CALCUL D'UNE BORNE POUR DP
!
    if (mode .eq. 0) then
        do n = 1, ndimsi
            v(n) = sigel(n)-mat(10)*a(n)
        end do
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        semax = sqrt(1.5d0*ddot(b_n, v, b_incx, v, b_incy))
        if (seeq .gt. semax) semax = seeq
        dpmax = (semax-d*mat(4))/1.5d0/(mat(2)+mat(11)*mat(10))
        goto 999
    end if
!
!
! -- CALCUL DE SIGMA ET DES CRITERES F ET G
!
    if (seeq .ne. 0.d0) then
        do n = 1, ndimsi
            s0(n) = se(n)/seeq
        end do
    else
        do n = 1, ndimsi
            s0(n) = 0.d0
        end do
        s0(4) = 1.d0/sqrt(1.5d0)
    end if
!
    if (mode .eq. 4) then
        do n = 1, ndimsi
            sig(n) = 1.5d0*mat(2)*dp*s0(n)
        end do
        goto 999
    end if
!
    do n = 1, ndimsi
        m(n) = b(n)+1.5d0*dp*s0(n)
        l(n) = a(n)+1.5d0*mat(11)*dp*s0(n)
    end do
!
    seq = seeq-1.5d0*(mat(2)+mat(11)*c)*dp
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    meq = sqrt(1.5d0*ddot(b_n, m, b_incx, m, b_incy))
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    leq = sqrt(1.5d0*ddot(b_n, l, b_incx, l, b_incy))
    if (meq .le. 0.d0) then
        q = mat(4)
    else
        q = mat(4)+mat(7)*meq**mat(5)
    end if
!
    f = seq-d*q
    if (dp .ne. 0.d0) f = f-mat(13)*dp**mat(12)*p**mat(14)
    g = c*leq+d*q-sp
    if (mode .eq. 1) goto 999
!
!
! -- DERIVEES DES CRITERES PAR RAPPORT A SIGMA PIC
!
    tmp = mat(8)*p/mat(11)*exp(-mat(8)*p*(1-sp/mat(11)))
    cds = mat(9)*tmp
    dds = -mat(6)*tmp
!
    do n = 1, ndimsi
        ads(n) = -epn(n)
        seds(n) = -cds*a(n)-c*ads(n)
    end do
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    seeqds = 1.5d0*ddot(b_n, se, b_incx, seds, b_incy)/seeq
    seqds = seeqds-1.5d0*mat(11)*cds*dp
!
    do n = 1, ndimsi
        s0ds(n) = (seds(n)*seq-seeqds*se(n))/seeq**2
        mds(n) = 1.5d0*dp*s0ds(n)
        lds(n) = ads(n)+1.5d0*mat(11)*dp*s0ds(n)
    end do
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    meqds = 1.5d0*ddot(b_n, m, b_incx, mds, b_incy)/meq
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    leqds = 1.5d0*ddot(b_n, l, b_incx, lds, b_incy)/leq
    qds = mat(7)*mat(5)*meq**(mat(5)-1)*meqds
!
    fds = seqds-dds*q-d*qds
    gds = cds*leq+c*leqds+dds*q+d*qds-1.d0
!
!
! -- DERIVEE DES CRITERES PAR RAPPORT A P
!
    if (mode .eq. 2 .or. mode .ge. 5) then
        tmp = -mat(8)*(1-sp/mat(11))*exp(-mat(8)*p*(1-sp/mat(11)))
        cdp = mat(9)*tmp
        ddp = -mat(6)*tmp
!
        do n = 1, ndimsi
            sedp(n) = -cdp*a(n)
        end do
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        seeqdp = 1.5d0*ddot(b_n, se, b_incx, sedp, b_incy)/seeq
        seqdp = seeqdp-1.5d0*(mat(2)+mat(11)*c+mat(11)*cdp*dp)
!
        do n = 1, ndimsi
            s0dp(n) = (sedp(n)*seeq-se(n)*seeqdp)/seeq**2
            mdp(n) = 1.5d0*(s0(n)+dp*s0dp(n))
            ldp(n) = mat(11)*mdp(n)
        end do
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        meqdp = 1.5d0*ddot(b_n, m, b_incx, mdp, b_incy)/meq
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        leqdp = 1.5d0*ddot(b_n, l, b_incx, ldp, b_incy)/leq
        qdp = mat(7)*mat(5)*meq**(mat(5)-1)*meqdp
!
        fdp = seqdp-ddp*q-d*qdp-mat(13)*(mat(14)*dp+mat(12)*p)*p**(mat(14)-1)*dp**(mat(12)-1)
        gdp = cdp*leq+c*leqdp+ddp*q+d*qdp
!
        if (mode .eq. 2) goto 999
    end if
!
!
! -- DERIVEES PAR RT A SIGEL ET CONSTRUCTION DE LA MATRICE TANGENTE
    if (mode .eq. 5 .or. mode .eq. 6) then
        tmp = 1.5d0/meq*1.5d0*dp/seeq
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pbs0 = 1.5d0*ddot(b_n, b, b_incx, s0, b_incy)
        do n = 1, ndimsi
            seqde(n) = 1.5d0*s0(n)
            meqde(n) = tmp*(b(n)-pbs0*s0(n))
            qde(n) = mat(7)*mat(5)*meq**(mat(5)-1)*meqde(n)
            fde(n) = seqde(n)-d*qde(n)
        end do
!
        if (mode .eq. 6) then
            tmp = 1.5d0/leq*1.5d0*mat(11)*dp/seeq
            b_n = to_blas_int(ndimsi)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            pas0 = 1.5d0*ddot(b_n, a, b_incx, s0, b_incy)
            do n = 1, ndimsi
                leqde(n) = tmp*(a(n)-pas0*s0(n))
                gde(n) = c*leqde(n)+d*qde(n)
            end do
        end if
!
        if (mode .eq. 5) then
            do n = 1, ndimsi
                pde(n) = -fde(n)/fdp
            end do
        else
            det = fdp*gds-gdp*fds
            do n = 1, ndimsi
                pde(n) = (fds*gde(n)-gds*fde(n))/det
            end do
        end if
!
        tmp = 1.5d0*dp/seeq
        do i = 1, ndimsi
            do j = 1, ndimsi
                tang(i, j) = 1.5d0*(s0(i)*pde(j)-tmp*s0(i)*s0(j))
            end do
        end do
        do n = 1, ndimsi
            tang(n, n) = tang(n, n)+tmp
        end do
        do i = 1, 3
            do j = 1, 3
                tang(i, j) = tang(i, j)-tmp/3.d0
            end do
        end do
!
        goto 999
    end if
!
!
! -- DERIVEES DES CRITERES PAR RAPPORT A XI
!
    if (mode .eq. 3) then
        do n = 1, ndimsi
            adx(n) = -sp*epndx(n)
            bdx(n) = -epndx(n)
            sedx(n) = -c*adx(n)
        end do
!
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        seeqdx = 1.5d0*ddot(b_n, se, b_incx, sedx, b_incy)/seeq
        seqdx = seeqdx
!
        do n = 1, ndimsi
            s0dx(n) = (sedx(n)*seq-seeqdx*se(n))/seeq**2
            mdx(n) = bdx(n)+1.5d0*dp*s0dx(n)
            ldx(n) = adx(n)+1.5d0*mat(11)*dp*s0dx(n)
        end do
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        meqdx = 1.5d0*ddot(b_n, m, b_incx, mdx, b_incy)/meq
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        leqdx = 1.5d0*ddot(b_n, l, b_incx, ldx, b_incy)/leq
        qdx = mat(7)*mat(5)*meq**(mat(5)-1)*meqdx
!
        fdx = seqdx-d*qdx
        gdx = c*leqdx+d*qdx
!
        goto 999
    end if
!
!
999 continue
end subroutine
