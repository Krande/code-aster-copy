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

subroutine te0226(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/vecma.h"
    character(len=16) :: option, nomte
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          COQUE 1D
!                          OPTION : 'MASS_MECA       '
!                          ELEMENT: MECXSE3
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
    character(len=8) :: elrefe, fami, poum
    integer(kind=8) :: icodre(1), kpg, spt
    real(kind=8) :: dfdx(3), r, rm, rf, rmf, poids, cour, nx, ny, h, vfi, vfj
    real(kind=8) :: matp(9, 9), matv(45), rho(1)
    integer(kind=8) :: nno, nnos, jgano, ndim, ipoids, ivf, idfdk, igeom, imate, icaco
    integer(kind=8) :: kp, npg, ii, jj, i, j, k, imatuu, kd1, kd2, kd3, ij1, ij2, ij3
    integer(kind=8) :: nddl, nvec, iacce, ivect
! ......................................................................
!
    call elref1(elrefe)
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
    nddl = 3*nno
    nvec = nddl*(nddl+1)/2
!
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCACOQU', 'L', icaco)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
    h = zr(icaco)
    rm = rho(1)*h
    rf = rho(1)*h**3/12.d0
!
    do k = 1, nvec
        matv(k) = 0.0d0
    end do
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                    cour, poids, nx, ny)
        if (nomte .eq. 'MECXSE3') then
            r = 0.0d0
            do i = 1, nno
                r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            end do
            poids = poids*r
            rmf = rf*(cour+nx/r)
        end if
!
        kd1 = 5
        kd2 = 3
        kd3 = 2
        do i = 1, 3*nno, 3
            kd1 = kd1+3*i-6
            kd2 = kd2+3*i-3
            kd3 = kd3+3*i
            ii = (i+2)/3
            do j = 1, i, 3
                jj = (j+2)/3
                ij1 = kd1+j-2
                ij2 = kd2+j-2
                ij3 = kd3+j-2
                vfi = zr(ivf+k+ii-1)
                vfj = zr(ivf+k+jj-1)
                matv(ij1) = matv(ij1)+vfi*vfj*poids*rm
                matv(ij2) = 0.0d0
                matv(ij2+1) = matv(ij1)
                matv(ij3) = matv(ij3)+vfi*vfj*poids*rmf*ny
                matv(ij3+1) = matv(ij3+1)-vfi*vfj*poids*rmf*nx
                matv(ij3+2) = matv(ij3+2)+vfi*vfj*poids*rf
            end do
!
            do j = 1, i-3, 3
                jj = (j+2)/3
                ij1 = kd1+j-2
                ij2 = kd2+j-2
                ij3 = kd3+j-2
                matv(ij1+1) = matv(ij2)
                matv(ij1+2) = matv(ij3)
                matv(ij2+2) = matv(ij3+1)
            end do
        end do
    end do
!
    if (option .eq. 'MASS_MECA') then
!
        call jevech('PMATUUR', 'E', imatuu)
!
        do i = 1, nvec
            zr(imatuu+i-1) = matv(i)
        end do
!
    else if (option .eq. 'M_GAMMA') then
!
        call jevech('PACCELR', 'L', iacce)
        call jevech('PVECTUR', 'E', ivect)
        call vecma(matv, nvec, matp, nddl)
        call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
!
    else
!C OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    end if
!
end subroutine
