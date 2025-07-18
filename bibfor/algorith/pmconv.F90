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
subroutine pmconv(r, rini, r1, inst, sigp, &
                  coef, iter, indimp, ds_conv, conver, &
                  itemax)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/pmimpr.h"
#include "asterfort/utmess.h"
#include "asterfort/GetResi.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    real(kind=8) :: r(12)
    real(kind=8) :: rini(12)
    real(kind=8) :: r1(12)
    real(kind=8) :: inst
    real(kind=8) :: sigp(6)
    real(kind=8) :: coef
    integer(kind=8) :: iter
    integer(kind=8) :: indimp(6)
    type(NL_DS_Conv), intent(in) :: ds_conv
    aster_logical :: conver
    aster_logical :: itemax
!
!-----------------------------------------------------------------------
!           OPERATEUR    CALC_POINT_MAT CALCUL D'ERREUR ET CONVERGENCE
!-----------------------------------------------------------------------
! IN   R      : RESIDU ACTUEL
! IN   RINI   : RESIDU INITIAL
! IN/OUT R1   : RESIDU PREMIERE ITERATION
! IN   INST   : INSTANT ACTUEL
! IN   SIGP   : CONTRAINTES ACTUELLES (POUR CONSTRUIRE LE DENOMINATEUR)
! IN   COEF   : COEF POUR ADIMENSIONNALISER LE PB
! IN   ITER   : NUMERO D'ITERATION
! In  ds_conv          : datastructure for convergence management
! OUT  ITEMAX : .TRUE. SI ITERATION MAXIMUM ATTEINTE
! OUT  CONVER : .TRUE. SI CONVERGENCE REALISEE
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ind, i, itmax, irela
    real(kind=8) :: resi_glob_rela, resi_glob_maxi
    real(kind=8) :: r8b(12)
    real(kind=8) :: ee, e1, e2, toler, e1ini, e2ini, er1, eini
    aster_logical :: l_rela
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
!
! - Get parameters
!
    call GetResi(ds_conv, type='RESI_GLOB_RELA', user_para_=resi_glob_rela, l_resi_test_=l_rela)
    call GetResi(ds_conv, type='RESI_GLOB_MAXI', user_para_=resi_glob_maxi)
!
!-----------------------------------------------------------------------
!     VERIFICATION DE LA CONVERGENCE EN DY  ET RE-INTEGRATION ?
!-----------------------------------------------------------------------
    conver = .false.
    e1 = 0.d0
    e2 = 0.d0
    e1ini = 0.d0
    e2ini = 0.d0
    er1 = 0.d0
    itemax = .false.
    if (iter .eq. 1) then
!        SAUVEGARDE DE R(DY0) POUR TEST DE CONVERGENCE
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigp, b_incx, r1(1), b_incy)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/coef, r1(1), b_incx)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, r(7), b_incx, r1(7), b_incy)
        do i = 1, 12
            er1 = max(er1, abs(r1(i)))
        end do
        if (er1 .le. r8prem()) then
            ee = er1
            ind = 4
            conver = .true.
            goto 999
        end if
    end if
!
    do i = 1, 6
        e1 = max(e1, abs(r(i)))
        e1ini = max(e1ini, abs(rini(i)))
        e1ini = max(e1ini, abs(r1(i)))
    end do
    do i = 7, 12
        e2 = max(e2, abs(r(i)))
        e2ini = max(e2ini, abs(rini(i)))
        e2ini = max(e2ini, abs(r1(i)))
    end do
    eini = max(e1ini, e2ini)
!
!     TEST RELATIF OU ABSOLU
    if (l_rela) then
        irela = 1
    else
        irela = 0
    end if
    if (irela .eq. 1) then
        toler = resi_glob_rela
        if (eini .gt. r8prem()) then
            e1 = e1/eini
            e2 = e2/eini
            ee = max(e1, e2)
            ind = 3
        end if
    else
        toler = resi_glob_maxi
        ee = max(e1, e2)
        ind = 4
    end if
    itemax = .false.
    itmax = ds_conv%iter_glob_maxi
!
    if (iter .lt. itmax) then
! -      NON CONVERGENCE ITERATION SUIVANTE
        if (ee .gt. toler) then
            conver = .false.
        else
            conver = .true.
        end if
    else
! -      NB ITERATION MAXIMUM ATTEINT SANS CONVERGENCE
        conver = .false.
        itemax = .true.
        call utmess('I', 'COMPOR2_5')
    end if
999 continue
!
    call pmimpr(ind, inst, indimp, r8b, iter, &
                r8b, r8b, r8b, 1, r8b, &
                ee, eini)
end subroutine
