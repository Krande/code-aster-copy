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

subroutine te0177(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/utmess.h"
#include "asterc/r8prem.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_ACOU'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    integer(kind=8) :: kp, i, j, k, ij, imattt, igeom, imate
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano
    real(kind=8) :: poids, r, cele_r, cele_i
    complex(kind=8) :: cele_c
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PMATTTC', 'E', imattt)
!
! - Get material properties for fluid
!
    call getFluidPara(zi(imate), cele_r_=cele_r, cele_i_=cele_i)

    cele_c = dcmplx(cele_r, cele_i)
!
! - Conditions on fluid parameters
!
    if ((abs(cele_r) .le. r8prem()) .and. (abs(cele_i) .le. r8prem())) then
        call utmess('F', 'FLUID1_7', sk='CELE_R + i*CELE_I')
    end if
!
! - Compute
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
        if (lteatt('AXIS', 'OUI')) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            end do
            poids = poids*r
        end if
!
        ij = imattt-1
        do i = 1, nno
            do j = 1, i
                ij = ij+1
                zc(ij) = zc(ij)+poids*((1.0d0, 0.0d0)/(cele_c**2))*zr(ivf+k+i-1)*zr(ivf+k+j-1)
            end do
        end do
    end do
!
end subroutine
