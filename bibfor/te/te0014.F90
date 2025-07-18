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

subroutine te0014(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
! Option: CHAR_MECA_ROTA_R
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: icodre(1)
    character(len=16) :: phenom
    real(kind=8) :: amm(81, 81), ft(81), x(27), y(27), z(27)
    real(kind=8) :: xi, xij
    real(kind=8) :: poids
    real(kind=8) :: rho(1), om1, om2, om3, omm, omo, rri
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: jgano, ndl, nno, kp, npg, ii, jj, i, j
    integer(kind=8) :: ndim, l, ic
    integer(kind=8) :: iret, nnos
    integer(kind=8) :: j_geom, j_rota, j_vect, j_mate, j_deplm, j_deplp
    real(kind=8) :: rota_speed, rota_axis(3), rota_cent(3)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'CHAR_MECA_ROTA_R')
!
! - Finite element parameters
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    ndl = 3*nno
    do i = 1, ndl
        do j = 1, ndl
            amm(i, j) = 0.d0
        end do
    end do
!
! - IN fields
!
    call jevech('PGEOMER', 'L', j_geom)
    call jevech('PMATERC', 'L', j_mate)
    call jevech('PROTATR', 'L', j_rota)
    call tecach('ONO', 'PDEPLMR', 'L', iret, iad=j_deplm)
    call tecach('ONO', 'PDEPLPR', 'L', iret, iad=j_deplp)
    rota_speed = zr(j_rota-1+1)
    rota_axis(1) = zr(j_rota-1+2)
    rota_axis(2) = zr(j_rota-1+3)
    rota_axis(3) = zr(j_rota-1+4)
    rota_cent(1) = zr(j_rota-1+5)
    rota_cent(2) = zr(j_rota-1+6)
    rota_cent(3) = zr(j_rota-1+7)
!
! - OUT fields
!
    call jevech('PVECTUR', 'E', j_vect)
!
! - Material
!
    call rccoma(zi(j_mate), 'ELAS', 1, phenom, icodre(1))
    call rcvalb('FPG1', 1, 1, '+', zi(j_mate), &
                ' ', phenom, 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre(1), 1)
!
! - Computation
!
    omm = rota_speed*rota_speed
    om1 = rota_speed*rota_axis(1)
    om2 = rota_speed*rota_axis(2)
    om3 = rota_speed*rota_axis(3)
    if (j_deplm .eq. 0 .or. j_deplp .eq. 0) then
        do i = 1, nno
            x(i) = zr(j_geom+3*i-3)-rota_cent(1)
            y(i) = zr(j_geom+3*i-2)-rota_cent(2)
            z(i) = zr(j_geom+3*i-1)-rota_cent(3)
        end do
    else
        do i = 1, nno
            x(i) = zr(j_geom+3*i-3)+zr(j_deplm+3*i-3)+zr(j_deplp+3*i-3)-rota_cent(1)
            y(i) = zr(j_geom+3*i-2)+zr(j_deplm+3*i-2)+zr(j_deplp+3*i-2)-rota_cent(2)
            z(i) = zr(j_geom+3*i-1)+zr(j_deplm+3*i-1)+zr(j_deplp+3*i-1)-rota_cent(3)
        end do
    end if
    do i = 1, nno
        omo = om1*x(i)+om2*y(i)+om3*z(i)
        ft(3*i-2) = omm*x(i)-omo*om1
        ft(3*i-1) = omm*y(i)-omo*om2
        ft(3*i) = omm*z(i)-omo*om3
    end do
!
! - Loop on point Gauss
!
    do kp = 1, npg
        l = (kp-1)*nno
        call dfdm3d(nno, kp, ipoids, idfde, zr(j_geom), &
                    poids)
        do i = 1, nno
            xi = rho(1)*poids*zr(ivf+l+i-1)
            ii = 3*(i-1)
            do j = 1, nno
                xij = xi*zr(ivf+l+j-1)
                jj = 3*(j-1)
                do ic = 1, 3
                    amm(ii+ic, jj+ic) = amm(ii+ic, jj+ic)+xij
                end do
            end do
        end do
    end do
!
    do i = 1, ndl
        rri = 0.d0
        do j = 1, ndl
            rri = rri+amm(i, j)*ft(j)
        end do
        amm(i, i) = rri
    end do
!
    do i = 1, ndl
        zr(j_vect+i-1) = amm(i, i)
    end do
!
end subroutine
