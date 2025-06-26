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
subroutine cvmjpl(mod, nmat, mater, timed, timef, &
                  epsd, deps, sigf, vinf, sigd, &
                  vind, nvi, nr, dsde)
! aslint: disable=W1306
    implicit none
!       VISCOCHABOCHE  :  MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT
!                         COHERENT A T OU T+DT
!       ----------------------------------------------------------------
!       IN  MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!           NR   :  DIMENSION DRDY
!           DRDY   :  MATRICE JACOBIENNE
!
!       DRDY  = ( DGDS  DGDX1  DGDX2  DGDP  DGDR  0   (DGDE3) )
!               ( DLDS  DLDX1  DLDX2  DLDP  DLDR  0   (DLDE3) )
!               ( DJDS  DJDX1  DJDX2  DJDP  DJDR  0   (DJDE3) )
!               ( DKDS  DKDX1  DKDX2  DKDP  DKDR  0   (DKDE3) )
!               ( DRDS  DRDX1  DRDX2  DRDP  DRDR  0   (DRDE3) )
!               ( 0     0      0      0     0     1   (0)     )
!               ((DQDS)(DQDX1)(DQDX2)(DQDP)(DQDR)(0)  (DQDE3) )
!                                                     ( SI IOPTIO = 1 )
!
!
!       DRDY  = ( DGDS  DGDX1  DGDX2  DGDP  DGDR  DGDQ  DGDXXI (DGDE3) )
!               ( DLDS  DLDX1  DLDX2  DLDP  DLDR  DLDQ  DLDXXI (DLDE3) )
!               ( DJDS  DJDX1  DJDX2  DJDP  DJDR  DJDQ  DJDXXI (DJDE3) )
!               ( DKDS  DKDX1  DKDX2  DKDP  DKDR  DKDQ  DKDXXI (DKDE3) )
!               ( DRDS  DRDX1  DRDX2  DRDP  DRDR  DRDQ  DRDXXI (DRDE3) )
!               ( DTDS  DTDX1  DTDX2  DTDP  DTDR  DTDQ  DTDXXI (DTDE3) )
!               ( DXIDS DXIDX1 DXIDX2 DXIDP DXIDR DXIDQ DXIDXI(DXIDE3))
!               ((DQDS)(DQDX1)(DQDX2)(DQDP)(DQDR)(DQDQ)(DQDXXI)(DQDE3) )
!                                                     ( SI IOPTIO = 2 )
!
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!       ----------------------------------------------------------------
#include "asterfort/cvmjac.h"
#include "asterfort/lcicma.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcprte.h"
#include "asterfort/mgauss.h"
#include "blas/daxpy.h"
    integer(kind=8) :: ndt, ndi, nmat, nr, nvi, iret
    integer(kind=8) :: ioptio, idnr, nopt
    real(kind=8) :: un, zero, det
    parameter(un=1.d0)
    parameter(zero=0.d0)
!
    real(kind=8) :: mater(nmat, 2)
!
    real(kind=8) :: dgds(6, 6), dgdx1(6, 6), dgdx2(6, 6), dgdr(6)
    real(kind=8) :: dlds(6, 6), dldx1(6, 6), dldx2(6, 6), dldr(6)
    real(kind=8) :: djds(6, 6), djdx1(6, 6), djdx2(6, 6), djdr(6)
    real(kind=8) :: dkds(6), dkdx1(6), dkdx2(6), dkdr
    real(kind=8) :: drds(6), drdx1(6), drdx2(6), drdr
    real(kind=8) :: dtds(6), dtdx1(6), dtdx2(6), dtdr
    real(kind=8) :: dxids(6, 6), dxidx1(6, 6), dxidx2(6, 6), dxidr(6)
    real(kind=8) :: dqds(6), dqdx1(6), dqdx2(6), dqdr
!
    real(kind=8) :: dgdq(6), dgdp(6), dgdxxi(6, 6), dgde3(6)
    real(kind=8) :: dldq(6), dldp(6), dldxxi(6, 6), dlde3(6)
    real(kind=8) :: djdq(6), djdp(6), djdxxi(6, 6), djde3(6)
    real(kind=8) :: dkdq, dkdp, dkdxxi(6), dkde3
    real(kind=8) :: drdq, drdp, drdxxi(6), drde3
    real(kind=8) :: dtdq, dtdp, dtdxxi(6), dtde3
    real(kind=8) :: dxidq(6), dxidp(6), dxidxi(6, 6), dxide3(6)
    real(kind=8) :: dqdq, dqdp, dqdxxi(6), dqde3
!
    real(kind=8) :: hookf(6, 6), dsde(6, 6), i6(6, 6)
!
    real(kind=8) :: matc(6, 6)
    real(kind=8) :: matd(6, 6), mate(6, 6), matf(6, 6)
    real(kind=8) :: mtmp(6, 6), mtmp1(6, 6), mtmp2(6, 6)
    real(kind=8) :: vtmp(6), vtmp1(6), vtmp2(6)
    real(kind=8) :: dkdset(6), dkdx1e(6), dkdx2e(6)
    real(kind=8) :: const1, const2, xx
!
    integer(kind=8) :: n1, n2, n3, n4, n5, n6, n7, n8, k
!
    character(len=8) :: mod
!
! DIMENSIONNEMENT DYNAMIQUE : TABLEAUX AUTOMATIQUES FORTRAN 90
    real(kind=8) :: drdy(nr, nr), yd(ndt+nvi), yf(ndt+nvi), dy(ndt+nvi)
    real(kind=8) :: timed, timef, sigf(*), vinf(*), sigd(*), vind(*)
    real(kind=8) :: epsd(*), deps(*)
    blas_int :: b_incx, b_incy, b_n
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/opti/ioptio, idnr
!       ----------------------------------------------------------------
    data i6/un, zero, zero, zero, zero, zero,&
     &                   zero, un, zero, zero, zero, zero,&
     &                   zero, zero, un, zero, zero, zero,&
     &                   zero, zero, zero, un, zero, zero,&
     &                   zero, zero, zero, zero, un, zero,&
     &                   zero, zero, zero, zero, zero, un/
!       ----------------------------------------------------------------
!
    nopt = 0
    do k = 1, 6
        dgde3(k) = 0.d0
    end do
    if (ioptio .eq. 2) nopt = idnr
!
    yd(1:ndt) = sigd(1:ndt)
    yf(1:ndt) = sigf(1:ndt)
    yd(ndt+1:ndt+nvi-1) = vind(1:nvi-1)
    yf(ndt+1:ndt+nvi-1) = vinf(1:nvi-1)
    dy(1:nr) = yf(1:nr)
    b_n = to_blas_int(nr)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, yd, b_incx, dy, &
               b_incy)
!
    call cvmjac(mod, nmat, mater, timed, timef, &
                yf, dy, nr, epsd, deps, &
                drdy)
!
!
! - RECUPERER LES SOUS-MATRICES BLOC
!
    n1 = 1
    n2 = ndt+1
    n3 = 2*ndt+1
    n4 = 3*ndt+1
    n5 = 3*ndt+2
    n6 = 3*ndt+3
    n7 = 3*ndt+4
    n8 = 3*ndt+4+nopt
!
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n1, n1, dgds, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n1, n2, dgdx1, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n1, n3, dgdx2, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n1, n4, dgdp, 6, 1, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n1, n5, dgdr, 6, 1, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n1, n6, dgdq, 6, 1, &
                1, 1)
!
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n2, n1, dlds, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n2, n2, dldx1, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n2, n3, dldx2, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n2, n4, dldp, 6, 1, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n2, n5, dldr, 6, 1, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n2, n6, dldq, 6, 1, &
                1, 1)
!
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n3, n1, djds, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n3, n2, djdx1, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, ndt, &
                n3, n3, djdx2, 6, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n3, n4, djdp, 6, 1, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n3, n5, djdr, 6, 1, &
                1, 1)
    call lcicma(drdy, nr, nr, ndt, 1, &
                n3, n6, djdq, 6, 1, &
                1, 1)
!
    call lcicma(drdy, nr, nr, 1, ndt, &
                n4, n1, dkds, 1, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, 1, ndt, &
                n4, n2, dkdx1, 1, 6, &
                1, 1)
!
    call lcicma(drdy, nr, nr, 1, ndt, &
                n4, n3, dkdx2, 1, 6, &
                1, 1)
!
    dkdp = drdy(n4, n4)
!
    dkdr = drdy(n4, n5)
!
    dkdq = drdy(n4, n6)
!
    dkdr = drdy(n4, n5)
!
    dkdq = drdy(n4, n6)
!
    call lcicma(drdy, nr, nr, 1, ndt, &
                n5, n1, drds, 1, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, 1, ndt, &
                n5, n2, drdx1, 1, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, 1, ndt, &
                n5, n3, drdx2, 1, 6, &
                1, 1)
!
    drdp = drdy(n5, n4)
!
    drdr = drdy(n5, n5)
!
    drdq = drdy(n5, n6)
!
    call lcicma(drdy, nr, nr, 1, ndt, &
                n6, n1, dtds, 1, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, 1, ndt, &
                n6, n2, dtdx1, 1, 6, &
                1, 1)
    call lcicma(drdy, nr, nr, 1, ndt, &
                n6, n3, dtdx2, 1, 6, &
                1, 1)
!
    dtdp = drdy(n6, n4)
!
    dtdr = drdy(n6, n5)
!
    dtdq = drdy(n6, n6)
!
!
    if (mod(1:6) .eq. 'C_PLAN') then
        call lcicma(drdy, nr, nr, ndt, 1, &
                    n1, n8, dgde3, 6, 1, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, 1, &
                    n2, n8, dlde3, 6, 1, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, 1, &
                    n3, n8, djde3, 6, 1, &
                    1, 1)
!
        dkde3 = drdy(n4, n8)
!
        drde3 = drdy(n5, n8)
!
        dtde3 = drdy(n6, n8)
!
!
        call lcicma(drdy, nr, nr, 1, ndt, &
                    n8, n1, dqds, 1, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, 1, ndt, &
                    n8, n2, dqdx1, 1, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, 1, ndt, &
                    n8, n3, dqdx2, 1, 6, &
                    1, 1)
!
        dqdp = drdy(n8, n4)
!
!
        dqdr = drdy(n8, n5)
!
!
        dqdq = drdy(n8, n6)
!
!
        dqde3 = drdy(n8, n8)
!
    end if
!
    if (ioptio .eq. 2) then
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n1, n7, dgdxxi, 6, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n2, n7, dldxxi, 6, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n3, n7, djdxxi, 6, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, 1, ndt, &
                    n4, n7, dkdxxi, 1, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, 1, ndt, &
                    n5, n7, drdxxi, 1, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, 1, ndt, &
                    n6, n7, dtdxxi, 1, 6, &
                    1, 1)
!
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n7, n1, dxids, 6, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n7, n2, dxidx1, 6, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n7, n3, dxidx2, 6, 6, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, 1, &
                    n7, n4, dxidp, 6, 1, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, 1, &
                    n7, n5, dxidr, 6, 1, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, 1, &
                    n7, n6, dxidq, 6, 1, &
                    1, 1)
        call lcicma(drdy, nr, nr, ndt, ndt, &
                    n7, n7, dxidxi, 6, 6, &
                    1, 1)
!
        if (mod(1:6) .eq. 'C_PLAN') then
!
            dqdq = drdy(n8, n6)
!
            call lcicma(drdy, nr, nr, 1, ndt, &
                        n8, n7, dqdxxi, 1, 6, &
                        1, 1)
            call lcicma(drdy, nr, nr, ndt, 1, &
                        n7, n8, dxide3, 6, 1, &
                        1, 1)
!
            dtde3 = drdy(n6, n8)
        end if
    end if
!
!       ----------------------------------------------------------------
!       L'OPTION 2 MODIFIE DKDS, DKDX1, DKDX2 CALCULES ICI SOUS LES NOMS
!       DE DKDSET, DKDX1E, ET DKDX2E
!       ----------------------------------------------------------------
    if (ioptio .eq. 2) then
        mtmp(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)
        call mgauss('NFVP', dxidxi, mtmp, 6, ndt, &
                    ndt, det, iret)
        vtmp2(1:ndt) = matmul(transpose(mtmp(1:ndt, 1:ndt)), dtdxxi(1:ndt))
        xx = dot_product(vtmp2(1:ndt), dxidp(1:ndt))
        const1 = dkdr/drdr*drdq
!
        if (const1 .ne. 0.d0) then
            const2 = const1/((dtdp-xx)*const1+dtdq*(dkdp-dkdr/drdr*drdp))
            const1 = 1.d0/const1
!
            vtmp(1:ndt) = matmul(transpose(dxids(1:ndt, 1:ndt)), vtmp2(1:ndt))
            vtmp1(1:ndt) = dtds(1:ndt)-vtmp(1:ndt)
            vtmp(1:ndt) = (const1*dtdq)*dkds(1:ndt)
            vtmp1(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
            dkdset(1:ndt) = const2*vtmp1(1:ndt)
!
            vtmp(1:ndt) = matmul(transpose(dxidx1(1:ndt, 1:ndt)), vtmp2(1:ndt))
            vtmp1(1:ndt) = dtdx1(1:ndt)-vtmp(1:ndt)
            vtmp(1:ndt) = (const1*dtdq)*dkdx1(1:ndt)
            vtmp1(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
            dkdx1e(1:ndt) = const2*vtmp1(1:ndt)
!
            vtmp(1:ndt) = matmul(transpose(dxidx2(1:ndt, 1:ndt)), vtmp2(1:ndt))
            vtmp1(1:ndt) = dtdx2(1:ndt)-vtmp(1:ndt)
            vtmp(1:ndt) = (const1*dtdq)*dkdx2(1:ndt)
            vtmp1(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
            dkdx2e(1:ndt) = const2*vtmp1(1:ndt)
!
        else
            const1 = 1.d0/(dkdp-dkdr/drdr*drdp)
            dkdset(1:ndt) = const1*dkds(1:ndt)
            dkdx1e(1:ndt) = const1*dkdx1(1:ndt)
            dkdx2e(1:ndt) = const1*dkdx2(1:ndt)
        end if
!
    else
        const1 = 1.d0/(dkdp-dkdr/drdr*drdp)
        dkdset(1:ndt) = const1*dkds(1:ndt)
        dkdx1e(1:ndt) = const1*dkdx1(1:ndt)
        dkdx2e(1:ndt) = const1*dkdx2(1:ndt)
    end if
!
! - E = ( DLDX2 - DLDP * DKDX2 ) * ( DJDX2 - DJDP * DKDX2 )-1
!
    call lcprte(djdp, dkdx2e, mtmp)
    mtmp(1:ndt, 1:ndt) = djdx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)
    call mgauss('NFVP', mtmp, mtmp1, 6, ndt, &
                ndt, det, iret)
    call lcprte(dldp, dkdx2e, mtmp)
    mtmp(1:ndt, 1:ndt) = dldx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mate(1:ndt, 1:ndt) = matmul(mtmp(1:ndt, 1:ndt), mtmp1(1:ndt, 1:ndt))
!
! - F = ( DJDX1 - DJDP * DKDX1 ) * ( DLDX1 - DLDP * DKDX1 )-1
!
    call lcprte(dldp, dkdx1e, mtmp)
    mtmp(1:ndt, 1:ndt) = dldx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)
    call mgauss('NFVP', mtmp, mtmp1, 6, ndt, &
                ndt, det, iret)
    call lcprte(djdp, dkdx1e, mtmp)
    mtmp(1:ndt, 1:ndt) = djdx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    matf(1:ndt, 1:ndt) = matmul(mtmp(1:ndt, 1:ndt), mtmp1(1:ndt, 1:ndt))
!
! - MATRICE C  TELLE QUE    DX1 = C * DSIG
!
    call lcprte(djdp, dkdx1e, mtmp)
    mtmp(1:ndt, 1:ndt) = djdx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = matmul(mate(1:ndt, 1:ndt), mtmp(1:ndt, 1:ndt))
    call lcprte(dldp, dkdx1e, mtmp)
    mtmp(1:ndt, 1:ndt) = dldx1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)-mtmp1(1:ndt, 1:ndt)
    mtmp2(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)
    call mgauss('NFVP', mtmp1, mtmp2, 6, ndt, &
                ndt, det, iret)
!
    call lcprte(djdp, dkdset, mtmp)
    mtmp(1:ndt, 1:ndt) = djds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = matmul(mate(1:ndt, 1:ndt), mtmp(1:ndt, 1:ndt))
    call lcprte(dldp, dkdset, mtmp)
    mtmp(1:ndt, 1:ndt) = dlds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = mtmp1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    matc(1:ndt, 1:ndt) = matmul(mtmp2(1:ndt, 1:ndt), mtmp(1:ndt, 1:ndt))
!
! - MATRICE D  TELLE QUE    DX2 = D * DSIG
!
    call lcprte(dldp, dkdx2e, mtmp)
    mtmp(1:ndt, 1:ndt) = dldx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = matmul(matf(1:ndt, 1:ndt), mtmp(1:ndt, 1:ndt))
    call lcprte(djdp, dkdx2e, mtmp)
    mtmp(1:ndt, 1:ndt) = djdx2(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)-mtmp1(1:ndt, 1:ndt)
    mtmp2(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)
    call mgauss('NFVP', mtmp1, mtmp2, 6, ndt, &
                ndt, det, iret)
!
    call lcprte(dldp, dkdset, mtmp)
    mtmp(1:ndt, 1:ndt) = dlds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = matmul(matf(1:ndt, 1:ndt), mtmp(1:ndt, 1:ndt))
    call lcprte(djdp, dkdset, mtmp)
    mtmp(1:ndt, 1:ndt) = djds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = mtmp1(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    matd(1:ndt, 1:ndt) = matmul(mtmp2(1:ndt, 1:ndt), mtmp(1:ndt, 1:ndt))
!
! - VTMP2 = DKDS + DKDX1 * C + DKDX2 * D
!
    vtmp2(1:ndt) = matmul(transpose(matd(1:ndt, 1:ndt)), dkdx2e(1:ndt))
    vtmp1(1:ndt) = matmul(transpose(matc(1:ndt, 1:ndt)), dkdx1e(1:ndt))
    vtmp2(1:ndt) = vtmp1(1:ndt)+vtmp2(1:ndt)
    vtmp2(1:ndt) = dkds(1:ndt)+vtmp2(1:ndt)
!
! - VTMP1 = DQDS + DQDX1 * C + DQDX2 * D - DQDP * VTMP2
!
    vtmp1(:) = 0.d0
    if (mod(1:6) .eq. 'C_PLAN') then
        vtmp(1:ndt) = dqdp*vtmp2(1:ndt)
        vtmp1(1:ndt) = dqds(1:ndt)-vtmp(1:ndt)
        vtmp(1:ndt) = matmul(transpose(matc(1:ndt, 1:ndt)), dqdx1(1:ndt))
        vtmp1(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
        vtmp2(1:ndt) = matmul(transpose(matd(1:ndt, 1:ndt)), dqdx2(1:ndt))
        vtmp1(1:ndt) = vtmp1(1:ndt)+vtmp(1:ndt)
    end if
!
! - MTMP  = DGDS + DGDX1 * C + DGDX2 * D - DGDP * VTMP2 - DGDE3 * VTMP1
!
    call lcprte(dgde3, vtmp1, mtmp1)
    call lcprte(dgdp, vtmp2, mtmp2)
    mtmp(1:ndt, 1:ndt) = mtmp1(1:ndt, 1:ndt)+mtmp2(1:ndt, 1:ndt)
    mtmp(1:ndt, 1:ndt) = dgds(1:ndt, 1:ndt)-mtmp(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = matmul(dgdx1(1:ndt, 1:ndt), matc(1:ndt, 1:ndt))
    mtmp(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
    mtmp1(1:ndt, 1:ndt) = matmul(dgdx2(1:ndt, 1:ndt), matd(1:ndt, 1:ndt))
    mtmp(1:ndt, 1:ndt) = mtmp(1:ndt, 1:ndt)+mtmp1(1:ndt, 1:ndt)
!
! - DSDE = (MTMP1)-1 * H
!
    mtmp1(1:ndt, 1:ndt) = i6(1:ndt, 1:ndt)
    call mgauss('NFVP', mtmp, mtmp1, 6, ndt, &
                ndt, det, iret)
    call lcopli('ISOTROPE', mod, mater(1, 1), hookf)
    dsde(1:ndt, 1:ndt) = matmul(mtmp1(1:ndt, 1:ndt), hookf(1:ndt, 1:ndt))
!
! - MATRICE DE COMPORTEMENT TANGENT:  SYMETRISATION DE DSDE
!
    mtmp(1:ndt, 1:ndt) = transpose(dsde(1:ndt, 1:ndt))
    dsde(1:ndt, 1:ndt) = dsde(1:ndt, 1:ndt)+mtmp(1:ndt, 1:ndt)
    dsde(1:ndt, 1:ndt) = 0.5d0*dsde(1:ndt, 1:ndt)
!
!
end subroutine
