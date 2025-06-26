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
subroutine nmhoff(ndim, imate, inst, epsm, deps, &
                  option, sigp, dsidep)
!
    implicit none
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvala.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dnrm2.h"
    integer(kind=8) :: ndim, imate
    real(kind=8) :: inst, deps(6), epsm(6), sigp(6), dsidep(6, 6)
    character(len=16) :: option
! -------------------------------------------------------------------
!     REALISE LA LOI DE NORTON-HOFF POUR LES ELEMENTS INCOMPRESSIBLES
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  INST    : INSTANT DU CALCUL
! IN  EPSM    : DEFORMATIONS A L'INSTANT PRECEDENT
! IN  DEPS    : INCREMENT DE DEFORMATIONS
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
!                                 RIGI_MECA_ELAS
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!
! --------------------------------------------------------------------
    aster_logical :: resi, rigi, elas, line
    integer(kind=8) :: k, l, ndimsi
    integer(kind=8) :: cod(1)
    real(kind=8) :: sy(1), m, am
    real(kind=8) :: eps(6), epsno
    real(kind=8) :: coef, rac23
    blas_int :: b_incx, b_incy, b_n
! -----------------------------------------------------------------
!
!
    call matfpe(-1)
!
! -- INITIALISATION
!
    ndimsi = 2*ndim
    rac23 = sqrt(2.d0/3.d0)
!
    resi = option .eq. 'RAPH_MECA' .or. option .eq. 'FULL_MECA'
    rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'FULL_MECA'
    elas = option .eq. 'RIGI_MECA_ELAS'
!
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, eps, b_incy)
    if (resi) then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, eps, &
                   b_incy)
    end if
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    epsno = dnrm2(b_n, eps, b_incx)
!
    line = inst .eq. 1 .or. epsno .eq. 0.d0
!
!
! -- CARACTERISTIQUES MATERIAU (SY) ET RIGIDITE
!
    call rcvala(imate, ' ', 'ECRO_LINE', 0, ' ', &
                [0.d0], 1, 'SY', sy, cod, &
                2)
    m = 1+10**(1-inst)
    am = sy(1)*rac23**m
    if (line) then
        coef = am
    else
        coef = am*epsno**(m-2)
    end if
!
!
! -- CALCUL DE SIGP
!
    if (resi) then
        do k = 1, ndimsi
            sigp(k) = coef*eps(k)
        end do
    end if
!
!
! -- CALCUL DE DSIDEP(6,6)
!
    if (rigi .or. elas) then
!
        call r8inir(36, 0.d0, dsidep, 1)
!
        do k = 1, ndimsi
            dsidep(k, k) = coef
        end do
!
        if (rigi .and. .not. line) then
            coef = coef*(m-2)/epsno**2
            do k = 1, ndimsi
                do l = 1, ndimsi
                    dsidep(k, l) = dsidep(k, l)+coef*eps(k)*eps(l)
                end do
            end do
        end if
    end if
!
    call matfpe(1)
!
end subroutine
