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
subroutine te0178(option, nomte)
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
    implicit none
!                          D'AMORTISSEMENT ACOUSTIQUE SUR DES ARETES
!                          D'ELEMENTS 2D
!                          OPTION : 'AMOR_ACOU'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/vff2dn.h"
#include "asterfort/getFluidPara.h"
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
!
    complex(kind=8) :: rhosz, cele_c
    character(len=16) :: option, nomte
    real(kind=8) :: poids, r, nx, ny, rho, alpha, q_alpha, cele_r, cele_i, onde_flui
    integer(kind=8) :: nno, kp, npg, ipoids, ivf, idfde, igeom, jv_amor
    integer(kind=8) :: imattt, i, j, ij, l, li, lj
    integer(kind=8) :: imate
    aster_logical :: laxi
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: jgano, mater, ndim, nnos
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PMATTTC', 'E', imattt)
    call jevech('PWATFLAC', 'L', jv_amor)
!
! - Get material properties for fluid
!
    mater = zi(imate)

    call getFluidPara(mater, rho_=rho, cele_r_=cele_r, cele_i_=cele_i, alpha_=alpha)

    cele_c = dcmplx(cele_r, cele_i)
!
! - Conditions on fluid parameters
!
    if ((1.d0-alpha) .ge. r8prem()) then
        q_alpha = (1.d0+alpha)/(1.d0-alpha)
    else
        call utmess('F', 'FLUID1_5', sk='ALPHA')
    end if

    if (rho .le. r8prem()) then
        call utmess('F', 'FLUID1_6', sk='RHO')
    end if

    if ((abs(cele_r) .le. r8prem()) .and. (abs(cele_i) .le. r8prem())) then
        call utmess('F', 'FLUID1_7', sk='CELE_R + i*CELE_I')
    end if

    if (zi(jv_amor-1+1) .eq. 1) then
        onde_flui = +1.d0
    else
        onde_flui = -1.d0
    end if
!
! - Output field
!

    rhosz = 1.d0/(q_alpha*cele_c)

!
! - Compute
!
    do kp = 1, npg
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
        if (laxi) then
            r = 0.d0
            do i = 1, nno
                l = (kp-1)*nno+i
                r = r+zr(igeom+2*i-2)*zr(ivf+l-1)
            end do
            poids = poids*r
        end if
        ij = imattt-1
        do i = 1, nno
            li = ivf+(kp-1)*nno+i-1
            do j = 1, i
                lj = ivf+(kp-1)*nno+j-1
                ij = ij+1
                zc(ij) = zc(ij)+poids*onde_flui*rhosz*zr(li)*zr(lj)
            end do
        end do
    end do
end subroutine
