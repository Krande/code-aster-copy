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

subroutine te0181(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/getFluidPara.h"
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DES MATRICES DE MASSE ELEMENTAIRE EN ACOUSTIQUE
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'MASS_ACOU '
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: jgano, nno, kp, npg, ij, i, j, imattt
    real(kind=8) :: poids, cele_r, cele_i
    complex(kind=8) :: cele_c
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: l, ndi, ndim, nnos
!-----------------------------------------------------------------------
    call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    ndi = nno*(nno+1)/2
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
    do i = 1, ndi
        zc(imattt-1+i) = (0.0d0, 0.0d0)
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg
        l = (kp-1)*nno
!
        call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
!
        do i = 1, nno
            do j = 1, i
                ij = (i-1)*i/2+j
                zc(imattt+ij-1) = zc(imattt+ij-1)+((1.0d0, 0.0d0)/(cele_c**2))*poids*zr(ivf&
                                  &+l+i-1)*zr(ivf+l+j-1)
            end do
        end do
    end do
!
end subroutine
