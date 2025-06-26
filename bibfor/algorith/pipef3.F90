! --------------------------------------------------------------------
! Copyright (C) 2007 NECS - BRUNO ZUBER   WWW.NECS.FR
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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
subroutine pipef3(ndim, nno, nddl, npg, lgpg, &
                  wref, vff, dfde, mate, geom, &
                  vim, ddepl, deplm, ddepl0, ddepl1, &
                  dtau, typmod, compor, copilo)
!
!
! aslint: disable=W1306
    implicit none
#include "asterc/r8vide.h"
#include "asterfort/nmfici.h"
#include "asterfort/pipeba.h"
#include "asterfort/pipetu.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    integer(kind=8) :: mate, npg, lgpg, nno, ndim, nddl
    real(kind=8) :: geom(nddl), vim(lgpg, npg), ddepl(nddl), deplm(nddl)
    real(kind=8) :: wref(npg), vff(nno, npg), dfde(2, nno, npg)
    real(kind=8) :: ddepl0(nddl), ddepl1(nddl), dtau
    character(len=8) :: typmod(2)
    character(len=16) :: compor
    real(kind=8) :: copilo(5, npg)
!
!-----------------------------------------------------------------------
!
!  PILOTAGE PRED_ELAS POUR LES ELEMENTS DE JOINT 3D
!
! IN  : GEOM, MATE, VIM, DDEPL, DEPLM, DDEPL0, DDELP1, DTAU, NPG
! OUT : COPILO
!-----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, kpg
    real(kind=8) :: up(nddl), ud(nddl), sup(3), sud(3), b(3, 60), poids
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
!
! DEPLACEMENT U(ETA) = UP + ETA * UD
!
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, deplm, b_incx, up, b_incy)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, ddepl, b_incx, up, &
               b_incy)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, ddepl0, b_incx, up, &
               b_incy)
    b_n = to_blas_int(nddl)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, ddepl1, b_incx, ud, b_incy)
!
! BOUCLE SUR LES POINTS DE GAUSS :
!
    do kpg = 1, npg
!
!      SAUT AU POINT DE GAUSS : SU(ETA) = SUP + ETA * SUD
        call nmfici(nno, nddl, wref(kpg), vff(1, kpg), dfde(1, 1, kpg), &
                    geom, poids, b)
!
        do i = 1, 3
            sup(i) = 0.d0
            sud(i) = 0.d0
            do j = 1, nddl
                sup(i) = sup(i)+b(i, j)*up(j)
                sud(i) = sud(i)+b(i, j)*ud(j)
            end do
        end do
!
!      INITIALISATION DES COEFFICIENTS DE PILOTAGE
        call r8inir(4, 0.d0, copilo(1, kpg), 1)
        copilo(5, kpg) = r8vide()
!
!      APPEL DU PILOTAGE PRED_ELAS SPECIFIQUE A LA LOI DE COMPORTEMENT
        if ((compor .eq. 'CZM_LIN_REG') .or. (compor .eq. 'CZM_EXP_REG')) then
            call pipeba(3, mate, sup, sud, vim(1, kpg), &
                        dtau, copilo(1, kpg))
        else if (compor .eq. 'CZM_TURON') then
            call pipetu(3, mate, sup, sud, vim(1, kpg), &
                        dtau, copilo(1, kpg))
        else
            call utmess('F', 'MECANONLINE_59')
        end if
!
    end do
!
end subroutine
