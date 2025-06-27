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
subroutine te0168(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/biline.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/vecma.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL MATRICE DE MASSE MECABLE
!                          OPTION : 'MASS_MECA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    integer(kind=8) :: icodre(1)
    real(kind=8) :: rho(1), coef, jacobi, en(3, 2)
    real(kind=8) :: matp(6, 6), matv(21), a
    integer(kind=8) :: nno, npg, k, kp, i, ii, jj, ki, ky, nddl, nvec, imatuu, lsect
    integer(kind=8) :: ipoids, ivf, iyty, igeom, imate, iacce, ivect
    integer(kind=8) :: ndim, nnos, jgano, idfdk
! ......................................................................
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
    call jevete('&INEL.CABPOU.YTY', 'L', iyty)
    nddl = 3*nno
    nvec = nddl*(nddl+1)/2
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
    call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
    call jevech('PCACABL', 'L', lsect)
    a = zr(lsect)
!
    k = 0
    do kp = 1, npg
        do i = 1, nno
            k = k+1
            en(i, kp) = zr(ivf-1+k)
        end do
    end do
!
    do k = 1, nvec
        matv(k) = 0.0d0
    end do
!
    do kp = 1, npg
        ky = (kp-1)*nddl*nddl
        jacobi = sqrt(biline(nddl, zr(igeom), zr(iyty+ky), zr(igeom)))
        coef = rho(1)*a*jacobi*zr(ipoids-1+kp)
        k = 0
        do ii = 1, nno
            do ki = 1, 3
                k = k+ki-3
                do jj = 1, ii
                    k = k+3
                    matv(k) = matv(k)+coef*en(ii, kp)*en(jj, kp)
                end do
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
!
        call vecma(matv, nvec, matp, nddl)
        call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
!
    else
!C OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    end if
!
end subroutine
