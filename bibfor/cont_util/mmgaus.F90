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

subroutine mmgaus(elem_type, type_inte, gauss_indx, xpg, ypg, &
                  gauss_weight_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/r8inir.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: elem_type
    integer(kind=8), intent(in) :: type_inte
    integer(kind=8), intent(in) :: gauss_indx
    real(kind=8), intent(out) :: xpg
    real(kind=8), intent(out) :: ypg
    real(kind=8), optional, intent(out) :: gauss_weight_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Continue method - Get integration points
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_type        : type of element
! In  type_inte        : type of integration scheme
! In  gauss_indx       : index of integration point
! Out xpg              : first parametric coordinates of integration point
! Out ypg              : second parametric coordinates of integration point
! Out gauss_weight     : weight of integration point
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: zgauss = 6
    integer(kind=8), parameter :: znpgse = 6
    integer(kind=8), parameter :: zseg = 2
    integer(kind=8), parameter :: znpgtr = 12
    integer(kind=8), parameter :: ztri = 3
    integer(kind=8), parameter :: zncots = 8
    integer(kind=8), parameter :: znpncs = 5
    integer(kind=8), parameter :: znpnct = 10
    real(kind=8) :: fpgseg(zgauss, znpgse, zseg)
    real(kind=8) :: fpgtri(zgauss, znpgtr, ztri)
    real(kind=8) :: pncseg(zncots, znpncs)
    real(kind=8) :: pnctri(zncots, znpnct)
    integer(kind=8) :: param
    integer(kind=8) :: i, j, h, n, incseg, jncseg
    real(kind=8) :: a, b, c, d, p1, p2, p3, gauss_weight
!
! --------------------------------------------------------------------------------------------------
!
    call r8inir(zgauss*znpgse*zseg, 0.d0, fpgseg, 1)
    call r8inir(zgauss*znpgtr*ztri, 0.d0, fpgtri, 1)
    call r8inir(zncots*znpncs, 0.d0, pncseg, 1)
    call r8inir(zncots*znpnct, 0.d0, pnctri, 1)
!
! - Gauss points for segment
!
    fpgseg(1, 1, 1) = 0.d0
    fpgseg(1, 1, 2) = 2.d0
    fpgseg(2, 1, 1) = +0.577350269189626d0
    fpgseg(2, 1, 2) = 1.d0
    fpgseg(2, 2, 1) = -fpgseg(2, 1, 1)
    fpgseg(2, 2, 2) = fpgseg(2, 1, 2)
!
    fpgseg(3, 1, 1) = -0.774596669241483d0
    fpgseg(3, 1, 2) = 0.555555555555556d0
    fpgseg(3, 2, 1) = 0.d0
    fpgseg(3, 2, 2) = 0.888888888888889d0
    fpgseg(3, 3, 1) = -fpgseg(3, 1, 1)
    fpgseg(3, 3, 2) = fpgseg(3, 1, 2)
!
    fpgseg(4, 1, 1) = +0.339981043584856d0
    fpgseg(4, 1, 2) = 0.652145154862546d0
    fpgseg(4, 2, 1) = -fpgseg(4, 1, 1)
    fpgseg(4, 2, 2) = fpgseg(4, 1, 2)
    fpgseg(4, 3, 1) = +0.861136311594053d0
    fpgseg(4, 3, 2) = 0.347854845137454d0
    fpgseg(4, 4, 1) = -fpgseg(4, 3, 1)
    fpgseg(4, 4, 2) = fpgseg(4, 3, 2)
!
    fpgseg(5, 1, 1) = -0.906179845938664d0
    fpgseg(5, 1, 2) = 0.236926885056189d0
    fpgseg(5, 2, 1) = -0.538469310105683d0
    fpgseg(5, 2, 2) = 0.478628670499366d0
    fpgseg(5, 3, 1) = 0.d0
    fpgseg(5, 3, 2) = 0.568888888888889d0
    fpgseg(5, 4, 1) = -fpgseg(5, 2, 1)
    fpgseg(5, 4, 2) = fpgseg(5, 2, 2)
    fpgseg(5, 5, 1) = -fpgseg(5, 1, 1)
    fpgseg(5, 5, 2) = fpgseg(5, 1, 2)
!
    fpgseg(6, 1, 1) = +0.238619186083197d0
    fpgseg(6, 1, 2) = 0.467913934572691d0
    fpgseg(6, 2, 1) = -fpgseg(6, 1, 1)
    fpgseg(6, 2, 2) = fpgseg(6, 1, 2)
    fpgseg(6, 3, 1) = +0.661209386466265d0
    fpgseg(6, 3, 2) = 0.360761573048139d0
    fpgseg(6, 4, 1) = -fpgseg(6, 3, 1)
    fpgseg(6, 4, 2) = fpgseg(6, 3, 2)
    fpgseg(6, 5, 1) = +0.932469514203152d0
    fpgseg(6, 5, 2) = 0.171324492379170d0
    fpgseg(6, 6, 1) = -fpgseg(6, 5, 1)
    fpgseg(6, 6, 2) = fpgseg(6, 5, 2)
!
! - Gauss points for triangle
!
    fpgtri(1, 1, 1) = 1.d0/3.d0
    fpgtri(1, 1, 2) = 1.d0/3.d0
    fpgtri(1, 1, 3) = 0.5d0
!
    a = 1.d0/6.d0
    b = 2.d0/3.d0
    p1 = 1.d0/6.d0
    fpgtri(2, 1, 1) = a
    fpgtri(2, 1, 2) = a
    fpgtri(2, 1, 3) = p1
    fpgtri(2, 2, 1) = b
    fpgtri(2, 2, 2) = a
    fpgtri(2, 2, 3) = p1
    fpgtri(2, 3, 1) = a
    fpgtri(2, 3, 2) = b
    fpgtri(2, 3, 3) = p1
!
    fpgtri(3, 1, 1) = 0.2d0
    fpgtri(3, 1, 2) = 0.2d0
    fpgtri(3, 1, 3) = 25.d0/96.d0
    fpgtri(3, 2, 1) = 0.6d0
    fpgtri(3, 2, 2) = 0.2d0
    fpgtri(3, 2, 3) = 25.d0/96.d0
    fpgtri(3, 3, 1) = 0.2d0
    fpgtri(3, 3, 2) = 0.6d0
    fpgtri(3, 3, 3) = 25.d0/96.d0
    fpgtri(3, 4, 1) = 1.d0/3.d0
    fpgtri(3, 4, 2) = 1.d0/3.d0
    fpgtri(3, 4, 3) = -27.d0/96.d0
!
    a = 0.445948490915965d0
    b = 0.091576213509771d0
    p1 = 0.111690794839005d0
    p2 = 0.054975871827661d0
    fpgtri(4, 1, 1) = b
    fpgtri(4, 1, 2) = b
    fpgtri(4, 1, 3) = p2
    fpgtri(4, 2, 1) = 1.d0-2.d0*b
    fpgtri(4, 2, 2) = b
    fpgtri(4, 2, 3) = p2
    fpgtri(4, 3, 1) = b
    fpgtri(4, 3, 2) = 1.d0-2.d0*b
    fpgtri(4, 3, 3) = p2
    fpgtri(4, 4, 1) = a
    fpgtri(4, 4, 2) = 1.d0-2.d0*a
    fpgtri(4, 4, 3) = p1
    fpgtri(4, 5, 1) = a
    fpgtri(4, 5, 2) = a
    fpgtri(4, 5, 3) = p1
    fpgtri(4, 6, 1) = 1.d0-2.d0*a
    fpgtri(4, 6, 2) = a
    fpgtri(4, 6, 3) = p1
!
    a = 0.470142064105115d0
    b = 0.101286507323456d0
    p1 = 0.066197076394253d0
    p2 = 0.062969590272413d0
    fpgtri(5, 1, 1) = 1.d0/3.d0
    fpgtri(5, 1, 2) = 1.d0/3.d0
    fpgtri(5, 1, 3) = 9.d0/80.d0
    fpgtri(5, 2, 1) = a
    fpgtri(5, 2, 2) = a
    fpgtri(5, 2, 3) = p1
    fpgtri(5, 3, 1) = 1.d0-2.d0*a
    fpgtri(5, 3, 2) = a
    fpgtri(5, 3, 3) = p1
    fpgtri(5, 4, 1) = a
    fpgtri(5, 4, 2) = 1.d0-2.d0*a
    fpgtri(5, 4, 3) = p1
    fpgtri(5, 5, 1) = b
    fpgtri(5, 5, 2) = b
    fpgtri(5, 5, 3) = p2
    fpgtri(5, 6, 1) = 1.d0-2.d0*b
    fpgtri(5, 6, 2) = b
    fpgtri(5, 6, 3) = p2
    fpgtri(5, 7, 1) = b
    fpgtri(5, 7, 2) = 1.d0-2.d0*b
    fpgtri(5, 7, 3) = p2
!
    a = 0.063089014491502d0
    b = 0.249286745170910d0
    c = 0.310352451033785d0
    d = 0.053145049844816d0
    p1 = 0.025422453185103d0
    p2 = 0.058393137863189d0
    p3 = 0.041425537809187d0
    fpgtri(6, 1, 1) = a
    fpgtri(6, 1, 2) = a
    fpgtri(6, 1, 3) = p1
    fpgtri(6, 2, 1) = 1.d0-2.d0*a
    fpgtri(6, 2, 2) = a
    fpgtri(6, 2, 3) = p1
    fpgtri(6, 3, 1) = a
    fpgtri(6, 3, 2) = 1.d0-2.d0*a
    fpgtri(6, 3, 3) = p1
    fpgtri(6, 4, 1) = b
    fpgtri(6, 4, 2) = b
    fpgtri(6, 4, 3) = p2
    fpgtri(6, 5, 1) = 1.d0-2.d0*b
    fpgtri(6, 5, 2) = b
    fpgtri(6, 5, 3) = p2
    fpgtri(6, 6, 1) = b
    fpgtri(6, 6, 2) = 1.d0-2.d0*b
    fpgtri(6, 6, 3) = p2
    fpgtri(6, 7, 1) = c
    fpgtri(6, 7, 2) = d
    fpgtri(6, 7, 3) = p3
    fpgtri(6, 8, 1) = d
    fpgtri(6, 8, 2) = c
    fpgtri(6, 8, 3) = p3
    fpgtri(6, 9, 1) = 1-c-d
    fpgtri(6, 9, 2) = c
    fpgtri(6, 9, 3) = p3
    fpgtri(6, 10, 1) = 1-c-d
    fpgtri(6, 10, 2) = d
    fpgtri(6, 10, 3) = p3
    fpgtri(6, 11, 1) = c
    fpgtri(6, 11, 2) = 1-c-d
    fpgtri(6, 11, 3) = p3
    fpgtri(6, 12, 1) = d
    fpgtri(6, 12, 2) = 1-c-d
    fpgtri(6, 12, 3) = p3
!
! - Newton-Cotes points for segment
!
    pncseg(3, 1) = 0.25d0
    pncseg(3, 2) = 0.75d0
!
    pncseg(4, 1) = 0.155555555555556d0
    pncseg(4, 2) = 0.711111111111111d0
    pncseg(4, 3) = 0.266666666666667d0
!
    pncseg(5, 1) = 0.131944444444444d0
    pncseg(5, 2) = 0.520833333333333d0
    pncseg(5, 3) = 0.347222222222222d0
!
    pncseg(6, 1) = 0.097619047619048d0
    pncseg(6, 2) = 0.514285714285714d0
    pncseg(6, 3) = 0.064285714285714d0
    pncseg(6, 4) = 0.647619047619048d0
!
    pncseg(7, 1) = 0.086921296296296d0
    pncseg(7, 2) = 0.414004629629630d0
    pncseg(7, 3) = 0.153125d0
    pncseg(7, 4) = 0.345949074074074d0
!
    pncseg(8, 1) = 0.069770723104057d0
    pncseg(8, 2) = 0.415379188712522d0
    pncseg(8, 3) = -0.065467372134039d0
    pncseg(8, 4) = 0.740458553791887d0
    pncseg(8, 5) = -0.320282186948854d0
!
! - Newton-cotes points for triangle
!
    pnctri(3, 1) = 0.016666666666667d0
    pnctri(3, 2) = 0.0375d0
    pnctri(3, 3) = 0.225d0
!
    pnctri(4, 1) = 0.d0
    pnctri(4, 2) = 0.044444444444445d0
    pnctri(4, 3) = -0.011111111111111d0
    pnctri(4, 4) = 0.088888888888889d0
!
    pnctri(5, 1) = 0.005456349206349d0
    pnctri(5, 2) = 0.012400793650794d0
    pnctri(5, 3) = 0.012400793650794d0
    pnctri(5, 4) = 0.099206349206349d0
    pnctri(5, 5) = 0.012400793650794d0
!
    pnctri(6, 1) = 0.d0
    pnctri(6, 2) = 0.021428571428572d0
    pnctri(6, 3) = -0.016071428571429d0
    pnctri(6, 4) = 0.042857142857143d0
    pnctri(6, 5) = 0.038095238095238d0
    pnctri(6, 6) = 0.042857142857143d0
    pnctri(6, 7) = -0.032142857142857d0
!
    pnctri(7, 1) = 0.002577160493827d0
    pnctri(7, 2) = 0.005765817901235d0
    pnctri(7, 3) = 0.006900077160494d0
    pnctri(7, 4) = 0.062195216049383d0
    pnctri(7, 5) = 0.005198688271605d0
    pnctri(7, 6) = -0.013233024691358d0
    pnctri(7, 7) = 0.086014660493827d0
    pnctri(7, 8) = 0.006616512345679d0
!
    pnctri(8, 1) = 0.d0
    pnctri(8, 2) = 0.012980599647266d0
    pnctri(8, 3) = -0.016507936507937d0
    pnctri(8, 4) = 0.024832451499118d0
    pnctri(8, 5) = 0.040070546737213d0
    pnctri(8, 6) = 0.029347442680776d0
    pnctri(8, 7) = -0.038201058201058d0
    pnctri(8, 8) = 0.023703703703704d0
    pnctri(8, 9) = -0.051075837742504d0
    pnctri(8, 10) = 0.051922398589065d0
!
!_______________________________________________________________________
!
! 'AUTO'
!
! TOUS LES SCHEMAS SONT DE TYPE TRAPEZE SAUF TR6/TR7/QU8/QU9
!
    if (type_inte .eq. 1) then
        if (elem_type(1:3) .eq. 'SE2') then
            if (gauss_indx .eq. 1) then
                xpg = -1.d0
                ypg = 0.d0
                gauss_weight = 1.d0
            else if (gauss_indx .eq. 2) then
                xpg = 1.d0
                ypg = 0.d0
                gauss_weight = 1.d0
            else
                ASSERT(.false.)
            end if
        else if (elem_type(1:3) .eq. 'SE3') then
            if (gauss_indx .eq. 1) then
                xpg = -1.d0
                ypg = 0.d0
                gauss_weight = 1.d0/3.d0
            else if (gauss_indx .eq. 2) then
                xpg = 1.d0
                ypg = 0.d0
                gauss_weight = 1.d0/3.d0
            else if (gauss_indx .eq. 3) then
                xpg = 0.d0
                ypg = 0.d0
                gauss_weight = 4.d0/3.d0
            else
                ASSERT(.false.)
            end if
        else if (elem_type(1:3) .eq. 'TR3') then
            if (gauss_indx .eq. 1) then
                xpg = 0.d0
                ypg = 0.d0
                gauss_weight = 1.d0/6.d0
            else if (gauss_indx .eq. 2) then
                xpg = 1.d0
                ypg = 0.d0
                gauss_weight = 1.d0/6.d0
            else if (gauss_indx .eq. 3) then
                xpg = 0.d0
                ypg = 1.d0
                gauss_weight = 1.d0/6.d0
            else
                ASSERT(.false.)
            end if
        else if ((elem_type(1:3) .eq. 'TR6') .or. (elem_type(1:3) .eq. 'TR7')) then
            ASSERT((gauss_indx .ge. 1) .and. (gauss_indx .le. 6))
            xpg = fpgtri(4, gauss_indx, 1)
            ypg = fpgtri(4, gauss_indx, 2)
            gauss_weight = fpgtri(4, gauss_indx, 3)
        else if ((elem_type(1:3) .eq. 'QU4')) then
            if (gauss_indx .eq. 1) then
                xpg = -1.d0
                ypg = -1.d0
                gauss_weight = 1.d0
            else if (gauss_indx .eq. 2) then
                xpg = 1.d0
                ypg = -1.d0
                gauss_weight = 1.d0
            else if (gauss_indx .eq. 3) then
                xpg = 1.d0
                ypg = 1.d0
                gauss_weight = 1.d0
            else if (gauss_indx .eq. 4) then
                xpg = -1.d0
                ypg = 1.d0
                gauss_weight = 1.d0
            else
                ASSERT(.false.)
            end if
        else if ((elem_type(1:3) .eq. 'QU8') .or. (elem_type(1:3) .eq. 'QU9')) then
            if (gauss_indx .eq. 1) then
                xpg = -1.d0
                ypg = -1.d0
                gauss_weight = 1.d0/9.d0
            else if (gauss_indx .eq. 2) then
                xpg = 1.d0
                ypg = -1.d0
                gauss_weight = 1.d0/9.d0
            else if (gauss_indx .eq. 3) then
                xpg = 1.d0
                ypg = 1.d0
                gauss_weight = 1.d0/9.d0
            else if (gauss_indx .eq. 4) then
                xpg = -1.d0
                ypg = 1.d0
                gauss_weight = 1.d0/9.d0
            else if (gauss_indx .eq. 5) then
                xpg = 0.d0
                ypg = -1.d0
                gauss_weight = 4.d0/9.d0
            else if (gauss_indx .eq. 6) then
                xpg = 1.d0
                ypg = 0.d0
                gauss_weight = 4.d0/9.d0
            else if (gauss_indx .eq. 7) then
                xpg = 0.d0
                ypg = 1.d0
                gauss_weight = 4.d0/9.d0
            else if (gauss_indx .eq. 8) then
                xpg = -1.d0
                ypg = 0.d0
                gauss_weight = 4.d0/9.d0
            else if (gauss_indx .eq. 9) then
                xpg = 0.d0
                ypg = 0.d0
                gauss_weight = 16.d0/9.d0
            else
                ASSERT(.false.)
            end if
        else
            ASSERT(.false.)
        end if
!_______________________________________________________________________
!
! 'GAUSS'
!
    else if (mod(type_inte, 10) .eq. 2) then
        param = type_inte/10
        if (elem_type(1:2) .eq. 'SE') then
            ASSERT((gauss_indx .ge. 1) .and. (gauss_indx .le. param))
            xpg = fpgseg(param, gauss_indx, 1)
            ypg = 0.d0
            gauss_weight = fpgseg(param, gauss_indx, 2)
        else if (elem_type(1:2) .eq. 'TR') then
            ASSERT((gauss_indx .ge. 1) .and. (gauss_indx .le. znpgtr))
            xpg = fpgtri(param, gauss_indx, 1)
            ypg = fpgtri(param, gauss_indx, 2)
            gauss_weight = fpgtri(param, gauss_indx, 3)
        else if (elem_type(1:2) .eq. 'QU') then
!           POINTS DE GAUSS ARRANGES EN LIGNE EN PARTANT DU BAS GAUCHE
            i = mod((gauss_indx-1), param)+1
            j = ((gauss_indx-1)/param)+1
!
            xpg = fpgseg(param, i, 1)
            ypg = fpgseg(param, j, 1)
            gauss_weight = fpgseg(param, i, 2)*fpgseg(param, j, 2)
        else
            ASSERT(.false.)
        end if
!_______________________________________________________________________
!
! 'SIMPSON'
!
    else if (mod(type_inte, 10) .eq. 3) then
        param = type_inte/10
!
! SEGMENTS
!
! EXEMPLE A 2 SUBDIVISIONS   NUMEROTATION            POIDS (X6)
!
!                         1---2---3---4---5      1---4---2---4---1
!
!   SUBDIVISIONS   | 1 | 2 | 3 | 4 |
! -----------------+---+---+---+---+
! NOMBRE DE POINTS | 3 | 5 | 7 | 9 |
!
        if (elem_type(1:2) .eq. 'SE') then
            n = gauss_indx-1
            xpg = -1.d0+n*(1.d0/param)
            ypg = 0.d0
!
            if ((n .eq. 0) .or. (n .eq. 2*param)) then
                gauss_weight = 1.d0/(param*3.d0)
            else
                if (mod(n, 2) .eq. 0) then
                    gauss_weight = 2.d0/(param*3.d0)
                else
                    gauss_weight = 4.d0/(param*3.d0)
                end if
            end if
!
! TRIANGLES
!
! EXEMPLE A 2 SUBDIVISIONS  NUMEROTATION            POIDS (X120)
!
!                           15                      1
!                           | \                     | \
!                           10  14                  4   4
!                           |     \                 |     \
!                           6   9   13              3---8---3
!                           |         \             | \     | \
!                           3   5   8   12          4   8   8   4
!                           |             \         |     \ |     \
!                           1---2---4---7--11       1---4---3---4---1
!
!    SUBDIVISIONS  | 1 | 2 | 3 | 4 |
! -----------------+---+---+---+---+
! NOMBRE DE POINTS | 6 | 15| 28| 45|
!
        else if (elem_type(1:2) .eq. 'TR') then
            h = 0
            n = gauss_indx
50          continue
            if (n .gt. 0) then
                h = h+1
                n = n-h
                goto 50
            end if
            h = h-1
            i = -n
            j = h+n
            h = 2*param
!
            xpg = 0.5d0*i/param
            ypg = 0.5d0*j/param
!
            if ((i .eq. 0) .or. (j .eq. 0) .or. (i+j .eq. h)) then
                if (i .eq. 0) then
                    if ((j .eq. 0) .or. (j .eq. h)) then
                        gauss_weight = 1.d0
                    else
                        if (mod(j, 2) .eq. 0) then
                            gauss_weight = 3.d0
                        else
                            gauss_weight = 4.d0
                        end if
                    end if
                else if (j .eq. 0) then
                    if (i .eq. h) then
                        gauss_weight = 1.d0
                    else
                        if (mod(i, 2) .eq. 0) then
                            gauss_weight = 3.d0
                        else
                            gauss_weight = 4.d0
                        end if
                    end if
                else
                    if (mod(j, 2) .eq. 0) then
                        gauss_weight = 3.d0
                    else
                        gauss_weight = 4.d0
                    end if
                end if
            else
                if ((mod(i, 2) .eq. 0) .and. (mod(j, 2) .eq. 0)) then
                    gauss_weight = 6.d0
                else
                    gauss_weight = 8.d0
                end if
            end if
            gauss_weight = gauss_weight/((param**2)*30.d0)
!
! QUADRANGLES
!
! EXEMPLE A 2 SUBDIVISIONS    NUMEROTATION           POIDS (X36)
!
!                          21--22--23--24--25     1---4---2---4---1
!                           |               |     |       |       |
!                          16  17  18  19  20     4  16   8  16   4
!                           |               |     |       |       |
!                          11  12  13  14  15     2---8---4---8---2
!                           |               |     |       |       |
!                           6   7   8   9  10     4  16   8  16   4
!                           |               |     |       |       |
!                           1---2---3---4---5     1---4---2---4---1
!
!   SUBDIVISIONS   | 1 | 2 | 3 | 4 |
! -----------------+---+---+---+---+
! NOMBRE DE POINTS | 9 | 25| 49| 81|
!
        else if (elem_type(1:2) .eq. 'QU') then
!
            i = mod((gauss_indx-1), (2*param+1))
            j = (gauss_indx-1)/(2*param+1)
!
            xpg = -1.d0+i*(1.d0/param)
            ypg = -1.d0+j*(1.d0/param)
!
            if ((i .eq. 0) .or. (j .eq. 0) .or. (i .eq. 2*param) .or. (j .eq. 2*param)) then
                if ((i .eq. 0) .or. (i .eq. 2*param)) then
                    if ((j .eq. 0) .or. (j .eq. 2*param)) then
                        gauss_weight = 1.d0/((param**2)*9.d0)
                    else
                        if (mod(j, 2) .eq. 0) then
                            gauss_weight = 2.d0/((param**2)*9.d0)
                        else
                            gauss_weight = 4.d0/((param**2)*9.d0)
                        end if
                    end if
                else
                    if (mod(i, 2) .eq. 0) then
                        gauss_weight = 2.d0/((param**2)*9.d0)
                    else
                        gauss_weight = 4.d0/((param**2)*9.d0)
                    end if
                end if
            else
                if ((mod(i, 2) .eq. 0) .and. (mod(j, 2) .eq. 0)) then
                    gauss_weight = 4.d0/((param**2)*9.d0)
                else if ((mod(i, 2) .eq. 1) .and. (mod(j, 2) .eq. 1)) then
                    gauss_weight = 16.d0/((param**2)*9.d0)
                else
                    gauss_weight = 8.d0/((param**2)*9.d0)
                end if
            end if
        else
            ASSERT(.false.)
        end if
!
!_______________________________________________________________________
!
!
! NCOTES
!
    else if (mod(type_inte, 10) .eq. 4) then
        param = type_inte/10
!
! SEGMENTS
!
! EXEMPLE D'ORDRE 4          NUMEROTATION            POIDS (X45)
!
!                         1---2---3---4---5      7---32--12--32--7
!
!      DEGRE       | 3 | 4 | 5 | 6 | 7 | 8 |
! -----------------+---+---+---+---+---+---+
! NOMBRE DE POINTS | 4 | 5 | 6 | 7 | 8 | 9 |
!
        if (elem_type(1:2) .eq. 'SE') then
            h = (param/2)+1
!
            if (gauss_indx .le. h) then
                incseg = gauss_indx
            else
                incseg = (param+1)-(gauss_indx-1)
            end if
!
            xpg = -1.d0+(gauss_indx-1)*(2.d0/param)
            ypg = 0.d0
            gauss_weight = pncseg(param, incseg)
!
! TRIANGLES
!
! EXEMPLE D'ORDRE 4   NUMEROTATION            POIDS (X90)
!
!                     15                      0
!                     | \                     | \
!                     10  14                  4   4
!                     |     \                 |     \
!                     6   9   13             -1   8  -1
!                     |         \             |         \
!                     3   5   8   12          4   8   8   4
!                     |             \         |             \
!                     1---2---4---7--11       0---4-(-1)--4---0
!
!      DEGRE       | 3 | 4 | 5 | 6 | 7 | 8 |
! -----------------+---+---+---+---+---+---+
! NOMBRE DE POINTS | 10| 15| 21| 28| 36| 45|
!
        else if (elem_type(1:2) .eq. 'TR') then
            h = 0
            n = gauss_indx
60          continue
            if (n .gt. 0) then
                h = h+1
                n = n-h
                goto 60
            end if
            i = -n
            j = h-1+n
!
            xpg = (i*1.d0)/param
            ypg = (j*1.d0)/param
!
            h = min(i, j, param-i-j)
            i = param-max(i, j, param-i-j)
!
            j = param/2
            if (i .le. j) then
                gauss_weight = pnctri(param, ((i+1)/2)*(i/2+1)+1+h)
            else
                gauss_weight = pnctri(param, ((i+1)/2)*(i/2+1)+1+h-(i-(param/2))* &
                                      (i-((param-1)/2)))
            end if
!
! QUADRANGLES
!
! EXEMPLE D'ORDRE 4
!    NUMEROTATION           POIDS (X2025)
!
!    21--22--23--24--25    49--224--84-224--49
!     |               |     |               |
!    16  17  18  19  20   224 1024 384 1024 224
!     |               |     |               |
!    11  12  13  14  15    84 384  144 384  84
!     |               |     |               |
!     6   7   8   9  10   224 1024 384 1024 224
!     |               |     |               |
!     1---2---3---4---5    49--224--84-224--49
!
!      DEGRE       | 3 | 4 | 5 | 6 | 7 | 8 |
! -----------------+---+---+---+---+---+---+
! NOMBRE DE POINTS | 16| 25| 36| 49| 64| 81|
!
        else if (elem_type(1:2) .eq. 'QU') then
            i = mod((gauss_indx-1), (param+1))+1
            j = (gauss_indx-1)/(param+1)+1
            h = (param/2)+1
!
            if (i .le. h) then
                incseg = i
            else
                incseg = (param+1)-(i-1)
            end if
            xpg = -1.d0+(i-1)*(2.d0/param)
!
            if (j .le. h) then
                jncseg = j
            else
                jncseg = (param+1)-(j-1)
            end if
            ypg = -1.d0+(j-1)*(2.d0/param)
!
            gauss_weight = pncseg(param, incseg)*pncseg(param, jncseg)
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
    end if
!
    if (present(gauss_weight_)) then
        gauss_weight_ = gauss_weight
    end if
!
end subroutine
