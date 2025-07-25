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

subroutine te0232(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
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
! Elements: COQUE_1D
! Option: CHAR_MECA_ROTA_R
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefe, fami, poum
    integer(kind=8) :: icodre(1), kpg, spt
    real(kind=8) :: zero, dfdx(3), nx, ny, poids, cour, rx, ry
    integer(kind=8) :: nno, kp, k, npg, i
    integer(kind=8) :: ipoids, ivf, idfdk
    integer(kind=8) :: jgano, ndim, nnos
    real(kind=8) :: rho(1)
    integer(kind=8) :: j_geom, j_rota, j_vect, j_mate, j_caco
    real(kind=8) :: rota_speed, rota_axis(3), rota_cent(3)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'CHAR_MECA_ROTA_R')
!
! - Finite element parameters
!
    call elref1(elrefe)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
! - IN fields
!
    call jevech('PGEOMER', 'L', j_geom)
    call jevech('PMATERC', 'L', j_mate)
    call jevech('PROTATR', 'L', j_rota)
    call jevech('PCACOQU', 'L', j_caco)
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
! - Checking
!

! AXE=Oy et CENTRE=ORIGINE
    if (abs(rota_axis(1)) .gt. r8miem() .or. abs(rota_axis(3)) .gt. r8miem()) then
        call utmess('F', 'CHARGES2_65')
    end if
    if (abs(rota_axis(2)) .le. r8miem()) then
        call utmess('F', 'CHARGES2_65')
    end if
    if (abs(rota_cent(1)) .gt. r8miem() .or. abs(rota_cent(2)) .gt. r8miem() .or. &
        abs(rota_cent(3)) .gt. r8miem()) then
        call utmess('F', 'CHARGES2_66')
    end if

!
! - Material
!
    zero = 0.d0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zi(j_mate), &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
!
! - Computation
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(j_geom), dfdx, &
                    cour, poids, nx, ny)
        poids = poids*rho(1)*rota_speed**2*zr(j_caco)
        rx = zero
        ry = zero
        do i = 1, nno
            rx = rx+zr(j_geom+2*i-2)*zr(ivf+k+i-1)
            ry = ry+zr(j_geom+2*i-1)*zr(ivf+k+i-1)
        end do
        poids = poids*rx
        do i = 1, nno
            zr(j_vect+3*i-3) = zr(j_vect+3*i-3)+poids*rota_axis(2)**2*rx*zr(ivf+k+i-1)
        end do
    end do
end subroutine
