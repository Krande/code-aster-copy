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

subroutine hujela(mod, mater, deps, sigd, sigf, iret)
    implicit none
!       INTEGRATION ELASTIQUE NON LINEAIRE DE LA LOI DE HUJEUX
!       IN  MOD    :  MODELISATION
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           SIGD   :  CONTRAINTE  A T
!           DEPS   :  INCREMENT DE DEFORMATION
!       OUT SIGF   :  CONTRAINTE A T+DT
!           IRET   :  CODE RETOUR DE  L'INTEGRATION DE LA LOI CJS
!                         IRET=0 => PAS DE PROBLEME
!                         IRET=1 => ECHEC
!       ---------------------------------------------------------------
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/hujci1.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, iret, i, j
    real(kind=8) :: coef, e, nu, al, demu, i1, n, pref
    real(kind=8) :: deps(6), dsig(6), sigd(6), sigf(6)
    real(kind=8) :: hook(6, 6), mater(22, 2)
    real(kind=8) :: zero, un, d13, deux, la, epsv, i1e
    real(kind=8) :: e1, e2, e3, nu12, nu13, nu23, g1, g2, g3, nu21, nu31, nu32
    real(kind=8) :: delta
    real(kind=8) :: piso, tole, c11, c12, c13, c22, c23, c33
    character(len=8) :: mod
    aster_logical :: tract
!
    common/tdim/ndt, ndi
!
    data zero/0.d0/
    data un/1.d0/
    data d13/0.33333333333334d0/
    data deux/2.d0/
    data tole/1.d-7/
!
!       ---------------------------------------------------------------
    pref = mater(8, 2)
    piso = 1.5d0*mater(21, 2)
    piso = zero
    n = mater(1, 2)
    i1e = d13*(sigd(1)+sigd(2)+sigd(3))
    epsv = deps(1)+deps(2)+deps(3)
    tract = .false.
!
    if (abs(n) .lt. r8prem()) then
!
        if (mater(17, 1) .eq. un) then
!
            i1 = i1e+d13*mater(1, 1)/(un-deux*mater(2, 1))*((i1e-piso)/pref)**n*epsv
!
        else if (mater(17, 1) .eq. deux) then
!
            e1 = mater(1, 1)
            e2 = mater(2, 1)
            e3 = mater(3, 1)
            nu12 = mater(4, 1)
            nu13 = mater(5, 1)
            nu23 = mater(6, 1)
            nu21 = mater(10, 1)
            nu31 = mater(11, 1)
            nu32 = mater(12, 1)
            delta = mater(13, 1)
!
            c11 = (un-nu23*nu32)*e1/delta
            c12 = (nu21+nu31*nu23)*e1/delta
            c13 = (nu31+nu21*nu32)*e1/delta
            c22 = (un-nu13*nu31)*e2/delta
            c23 = (nu32+nu31*nu12)*e2/delta
            c33 = (un-nu21*nu12)*e3/delta
!
            i1 = (c11+c12+c13)*deps(1)+(c12+c22+c23)*deps(2)+(c13+c23+c33)*deps(3)
            i1 = i1e+d13*i1
!
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if ((i1-piso)/pref .lt. tole) then
            tract = .true.
            goto 5
        end if
        goto 30
!
    end if
!
!
!--->  CALCUL DE I1=TR(SIG) A T+DT PAR METHODE DE LA SECANTE
!      OU EXPLICITEMENT SI NIVEAU HUJEUX
    call hujci1(mater, deps, sigd, i1, tract, iret)
!
    if (iret .eq. 1) goto 999
!
    if (mater(17, 1) .eq. un) then
!
        i1e = i1e+d13*mater(1, 1)/(un-deux*mater(2, 1))*epsv*((i1e-piso)/pref)**n
!
    else if (mater(17, 1) .eq. deux) then
!
        e1 = mater(1, 1)*((i1e-piso)/pref)**n
        e2 = mater(2, 1)*((i1e-piso)/pref)**n
        e3 = mater(3, 1)*((i1e-piso)/pref)**n
        nu12 = mater(4, 1)
        nu13 = mater(5, 1)
        nu23 = mater(6, 1)
        nu21 = mater(13, 1)
        nu31 = mater(14, 1)
        nu32 = mater(15, 1)
        delta = mater(16, 1)
!
        c11 = (un-nu23*nu32)*e1/delta
        c12 = (nu21+nu31*nu23)*e1/delta
        c13 = (nu31+nu21*nu32)*e1/delta
        c22 = (un-nu13*nu31)*e2/delta
        c23 = (nu32+nu31*nu12)*e2/delta
        c33 = (un-nu21*nu12)*e3/delta
!
        coef = (c11+c12+c13)*deps(1)+(c12+c22+c23)*deps(2)+(c13+ &
                                                            c23+c33)*deps(3)
        i1e = i1e+d13*coef
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if ((i1-piso)/pref .lt. tole) tract = .true.
!
!--->  EN CAS D'ENTREE EN TRACTION, LES CONTRAINTES SONT
!      RAMENEES SUR L'AXE HYDROSTATIQUE A DES VALEURS FAIBLES
!      ( JUSTE AU-DELA DE LA TOLERANCE )
5   continue
    if (tract) then
        do i = 1, ndi
            sigf(i) = pref*tole*1.01d0+piso
        end do
        do i = ndi+1, ndt
            sigf(i) = zero
        end do
        goto 999
    end if
30  continue
!
!
!---> CALCUL DU COEF  (-----------)**N ET MODULE_YOUNG A T+DT
    hook(:, :) = zero
!
    coef = ((i1-piso)/pref)**n
!
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
!
        if (mater(17, 1) .eq. un) then
!
            e = mater(1, 1)*coef
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
            e1 = mater(1, 1)*coef
            e2 = mater(2, 1)*coef
            e3 = mater(3, 1)*coef
            nu12 = mater(4, 1)
            nu13 = mater(5, 1)
            nu23 = mater(6, 1)
            g1 = mater(7, 1)*coef
            g2 = mater(8, 1)*coef
            g3 = mater(9, 1)*coef
            nu21 = mater(13, 1)
            nu31 = mater(14, 1)
            nu32 = mater(15, 1)
            delta = mater(16, 1)
!
            hook(1, 1) = (un-nu23*nu32)*e1/delta
            hook(1, 2) = (nu21+nu31*nu23)*e1/delta
            hook(1, 3) = (nu31+nu21*nu32)*e1/delta
            hook(2, 2) = (un-nu13*nu31)*e2/delta
            hook(2, 3) = (nu32+nu31*nu12)*e2/delta
            hook(3, 3) = (un-nu21*nu12)*e3/delta
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
!
    else if (mod(1:6) .eq. 'C_PLAN' .or. mod(1:2) .eq. '1D') then
!
        call utmess('F', 'COMPOR1_4')
!
    end if
!
!
!--->   INCREMENTATION DES CONTRAINTES  SIGF = SIGD + HOOK DEPS
    sigf(:) = 0.d0
    dsig(1:ndt) = matmul(hook(1:ndt, 1:ndt), deps(1:ndt))
    sigf(1:ndt) = sigd(1:ndt)+dsig(1:ndt)
!
999 continue
end subroutine
