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
subroutine lcmmjp(mod, nmat, mater, timed, timef, &
                  comp, nbcomm, cpmono, pgl, nfs, &
                  nsg, toutms, hsr, nr, nvi, &
                  sigd, itmax, toler, vinf, vind, &
                  dsde, drdy, option, iret)
! aslint: disable=W1306,W1504
    implicit none
!     COMPORTEMENT MONOCRISTALLIN
!                :  MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT
!                   COHERENT A T+DT? en hpp et gdef
!     ----------------------------------------------------------------
!     IN  MOD    :  TYPE DE MODELISATION
!         NMAT   :  DIMENSION MATER
!         MATER  :  COEFFICIENTS MATERIAU
!         TIMED  :  ISTANT PRECEDENT
!         TIMEF  :  INSTANT ACTUEL
!         COMP   :  NOM COMPORTEMENT
!         NBCOMM :  INCIDES DES COEF MATERIAU
!         CPMONO :  NOM DES COMPORTEMENTS
!         PGL    :  MATRICE DE PASSAGE
!         TOUTMS :  TENSEURS D'ORIENTATION
!         HSR    :  MATRICE D'INTERACTION
!         NVI    :  NOMBRE DE VARIABLES INTERNES
!         NR     :  DIMENSION DU SYSTEME A RESOUDRE
!         ITMAX  :  ITER_INTE_MAXI
!         TOLER  :  RESI_INTE
!         VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT T
!         VINF   :  VARIABLES INTERNES A T+DT
!         DRDY   :  MATRICE JACOBIENNE
!         OPTION :  OPTION DE CALCUL MATRICE TANGENTE
!     OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!                   DSDE = INVERSE(Y0-Y1*INVERSE(Y3)*Y2)
!         IRET   :  CODE RETOUR
!     ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/lcmmja.h"
#include "asterfort/lcmmkg.h"
#include "asterfort/mgauss.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/lcafyd.h"
#include "blas/dcopy.h"
    integer(kind=8) :: ndt, ndi, nmat, nvi, itmax, nfs, nsg
    integer(kind=8) :: k, j, nr, iret, ns, nbcomm(nmat, 3)
! DIMENSIONNEMENT DYNAMIQUE
    real(kind=8) :: drdy(nr, nr), dsde(6, *), kyl(6, 6), det, i6(6, 6)
    real(kind=8) :: zinv(6, 6)
    real(kind=8) :: toler, mater(*), yf(nr), dy(nr), un, zero, timed, timef
    real(kind=8) :: pgl(3, 3), sigd(6)
    real(kind=8) :: z0(6, 6), z1(6, (nr-ndt))
    real(kind=8) :: z2((nr-ndt), 6), z3((nr-ndt), (nr-ndt))
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg)
    real(kind=8) :: vind(*), vinf(*), df(9), yd(nr)
    character(len=8) :: mod
    character(len=16) :: comp(*), option
    character(len=24) :: cpmono(5*nmat+1)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    common/tdim/ndt, ndi
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
    data i6/un, zero, zero, zero, zero, zero,&
     &                 zero, un, zero, zero, zero, zero,&
     &                 zero, zero, un, zero, zero, zero,&
     &                 zero, zero, zero, un, zero, zero,&
     &                 zero, zero, zero, zero, un, zero,&
     &                 zero, zero, zero, zero, zero, un/
!
! -  INITIALISATION
!
    ns = nr-ndt
    iret = 0
!
!     RECALCUL DE LA DERNIERE MATRICE JACOBIENNE
    if (option .eq. 'RIGI_MECA_TANG') then
!
        call r8inir(nr, 0.d0, dy, 1)
        call r8inir(9, 0.d0, df, 1)
!       call r8inir(ndt, 0.d0, sigd, 1)
        call lcafyd(comp, mater, mater, nbcomm, cpmono, &
                    nmat, mod, nvi, vind, &
                    sigd, nr, yd)
        b_n = to_blas_int(nr)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, yd, b_incx, yf, b_incy)
        call lcmmja(mod, nmat, mater, timed, timef, &
                    itmax, toler, nbcomm, cpmono, pgl, &
                    nfs, nsg, toutms, hsr, nr, &
                    nvi, vind, df, yf, yd, &
                    dy, drdy, iret)
        if (iret .gt. 0) goto 999
    end if
!
! - RECUPERER LES SOUS-MATRICES BLOC
!
    do k = 1, 6
        do j = 1, 6
            z0(k, j) = drdy(k, j)
        end do
    end do
    do k = 1, 6
        do j = 1, ns
            z1(k, j) = drdy(k, ndt+j)
        end do
    end do
    do k = 1, ns
        do j = 1, 6
            z2(k, j) = drdy(ndt+k, j)
        end do
    end do
    do k = 1, ns
        do j = 1, ns
            z3(k, j) = drdy(ndt+k, ndt+j)
        end do
    end do
!     Z2=INVERSE(Z3)*Z2
!     CALL MGAUSS ('NCSP',Z3, Z2, NS, NS, 6, DET, IRET )
    call mgauss('NCWP', z3, z2, ns, ns, &
                6, det, iret)
    if (iret .gt. 0) goto 999
!
!     KYL=Z1*INVERSE(Z3)*Z2
    call promat(z1, 6, 6, ns, z2, &
                ns, ns, 6, kyl)
!
!     Z0=Z0+Z1*INVERSE(Z3)*Z2
    do k = 1, 6
        do j = 1, 6
            z0(k, j) = z0(k, j)-kyl(k, j)
        end do
    end do
!
    b_n = to_blas_int(36)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, i6, b_incx, zinv, b_incy)
!     CALL MGAUSS ('NCSP',Z0, ZINV, 6, 6, 6, DET, IRET )
    call mgauss('NCWP', z0, zinv, 6, 6, &
                6, det, iret)
    if (iret .gt. 0) goto 999
!
    if (gdef .eq. 0) then
!
!        DSDE = INVERSE(Z0-Z1*INVERSE(Z3)*Z2)
!
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zinv, b_incx, dsde, b_incy)
!
    else
!
        call lcmmkg(zinv, nvi, vind, vinf, nmat, &
                    mater, mod, nr, dsde)
!
    end if
999 continue
end subroutine
