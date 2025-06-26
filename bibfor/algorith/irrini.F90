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
subroutine irrini(fami, kpg, ksp, typess, essai, &
                  mod, nmat, materf, yd, deps, &
                  dy)
!
    implicit none
!
#include "asterfort/lcdevi.h"
#include "asterfort/lcopli.h"
#include "asterfort/rcvarc.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: typess, nmat, kpg, ksp
    real(kind=8) :: essai, materf(nmat, 2), yd(*), deps(6), dy(*)
    character(len=8) :: mod
    character(len=*) :: fami
!
! person_in_charge: jean-luc.flejou at edf.fr
!       IRRAD3M : CALCUL SOLUTION ESSAI DY = ( DSIG DX1 DX2 DP (DEPS3))
!                               AVEC     Y  = ( SIG  X1  X2  P  (EPS3))
!       IN  ESSAI  :  VALEUR DE LA SOLUTION D ESSAI
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           YD     :  VARIABLES A T   = ( SIG  VIN  (EPS3)  )
!       VAR DEPS   :  INCREMENT DE DEFORMATION
!           TYPESS :  TYPE DE SOLUTION D ESSAI
!                               0 = NUL(0)
!                               1 = ELASTIQUE
!                               2 = EXPLICITE (=-1 INITIALEMENT)
!                               3 = ESSAI
!       OUT DY     :  SOLUTION ESSAI  = ( DSIG DVIN (DEPS3) )
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
    real(kind=8) :: hook(6, 6), dev(6), s, dfds(6), vtmp1(6), vtmp2(6), dsig(6)
    real(kind=8) :: dphi, id3d(6), nun, sig(6), p, etai
    real(kind=8) :: k, n, p0, ai0, etais, ag, irrad, irraf, zetaf, zetag
    real(kind=8) :: detai, dpi, dp, dg, yy, xx, zz
    real(kind=8) :: penpe, pe, pk
    integer(kind=8) :: ndt, ndi, iret, i
    blas_int :: b_incx, b_incy, b_n
    data id3d/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
    if (typess .eq. -1) typess = 2
    sig(1:ndt) = yd(1:ndt)
    p = yd(ndt+1)
    etai = yd(ndt+2)
!
!     PARAMETRES MATERIAUX
    ai0 = materf(4, 2)
    etais = materf(5, 2)
    ag = materf(6, 2)
    k = materf(7, 2)
    n = materf(8, 2)
    p0 = materf(9, 2)
    zetaf = materf(12, 2)
    penpe = materf(13, 2)
    pk = materf(14, 2)
    pe = materf(15, 2)
    zetag = materf(17, 2)
!
!     POUR LES CONTRAINTES PLANES
    nun = materf(2, 1)/(1.d0-materf(2, 1))
!
    typess = 1
!     SOLUTION NULLE ( TYPESS=0) OU ELASTIQUE ( TYPESS=1)
    if (typess .eq. 0 .or. typess .eq. 1) then
        dy(1:ndt+4) = 0.d0
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = 0.d0
            dy(ndt+5) = 0.d0
        end if
!
        if (typess .eq. 1) then
            call lcopli('ISOTROPE', mod, materf(1, 1), hook)
            dy(1:ndt) = matmul(hook(1:ndt, 1:ndt), deps(1:ndt))
        end if
!
!     SOLUTION EXPLICITE
    else if (typess .eq. 2) then
        call lcopli('ISOTROPE', mod, materf(1, 1), hook)
        call rcvarc('F', 'IRRA', '-', fami, kpg, &
                    ksp, irrad, iret)
        call rcvarc('F', 'IRRA', '+', fami, kpg, &
                    ksp, irraf, iret)
!        ARRET DANS IRRMAT SI  IRRAD .GT. IRRAF*1.00001
        if (irrad .gt. irraf) then
            dphi = 0.0d0
        else
            dphi = irraf-irrad
        end if
!
        call lcdevi(sig, dev)
        b_n = to_blas_int(ndt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        s = ddot(b_n, dev, b_incx, dev, b_incy)
        s = sqrt(1.5d0*s)
!
!        DETAI
        detai = zetaf*s*dphi
!        DPI
        if ((etai+detai) .lt. etais) then
            dpi = 0.d0
        else if (etai .ge. etais) then
            dpi = ai0*detai
        else
            dpi = ai0*(detai-etais+etai)
        end if
!        DG
        dg = ag*dphi*zetag
!        DP
        if (s .eq. 0.d0) then
            dp = 0.d0
            do i = 1, 6
                dfds(i) = 0.d0
            end do
        else
            dfds(1:ndt) = (1.5d0/s)*dev(1:ndt)
            vtmp1(1:ndt) = dpi*dfds(1:ndt)
            vtmp1(1:ndt) = deps(1:ndt)-vtmp1(1:ndt)
            vtmp2(1:ndt) = dg*id3d(1:ndt)
            vtmp1(1:ndt) = vtmp1(1:ndt)-vtmp2(1:ndt)
            vtmp1(1:ndt) = matmul(hook(1:ndt, 1:ndt), vtmp1(1:ndt))
            yy = dot_product(dfds(1:ndt), vtmp1(1:ndt))
!
            if (p .lt. pk) then
                zz = 0.0d0
            else if (p .lt. pe) then
                zz = penpe
            else
                zz = n*k*(p+p0)**(n-1.0d0)
            end if
            vtmp1(1:ndt) = matmul(hook(1:ndt, 1:ndt), dfds(1:ndt))
            xx = dot_product(dfds(1:ndt), vtmp1(1:ndt))
!
            xx = xx+zz
!
            dp = yy/xx
        end if
!
!        (DEPS(3))
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = nun*((dp+dpi)*(dfds(1)+dfds(2))+2.d0*dg-deps(1)-deps(2))+dfds(3)*(dp+dpi &
                                                                                        )+dg
        end if
!
!        DSIG
        vtmp1(1:ndt) = ((dpi+dp))*dfds(1:ndt)
        vtmp1(1:ndt) = deps(1:ndt)-vtmp1(1:ndt)
        vtmp2(1:ndt) = dg*id3d(1:ndt)
        vtmp1(1:ndt) = vtmp1(1:ndt)-vtmp2(1:ndt)
        dsig(1:ndt) = matmul(hook(1:ndt, 1:ndt), vtmp1(1:ndt))
!        DY
        dy(1:ndt) = dsig(1:ndt)
        dy(ndt+1) = dp
        dy(ndt+2) = detai
        dy(ndt+3) = dpi
        dy(ndt+4) = dg
        if (mod(1:6) .eq. 'C_PLAN') then
            dy(ndt+5) = deps(3)
            dy(3) = 0.d0
        end if
!
! - SOLUTION INITIALE = VALEUR ESSAI POUR TOUTES LES COMPOSANTES
!
    else if (typess .eq. 3) then
        dy(1:ndt+4) = essai
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = essai
            dy(3) = 0.d0
        end if
    end if
!
end subroutine
