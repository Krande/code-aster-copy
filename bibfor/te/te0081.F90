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

subroutine te0081(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bmatmc.h"
#include "asterfort/btdbmc.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/get_elas_id.h"
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D
! Option: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, igau, igeom, imate, imatuu, j, k
    integer(kind=8) :: nbinco, nbsig
    real(kind=8) :: b(486), btdb(81, 81), d(16), jacgau
    real(kind=8) :: angl_naut(3), instan, nharm
    integer(kind=8) :: ndim, nno, nnos, npg1, ipoids, ivf, idfde
    character(len=4) :: fami
    integer(kind=8) :: elas_id
!
! --------------------------------------------------------------------------------------------------
!
!
! - Finite element informations
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Initializations
!
    ndim = 2
    instan = r8vide()
    nbinco = ndim*nno
    nharm = 0.d0
    btdb(:, :) = 0.d0
!
! - Number of stress components
!
    nbsig = nbsigm()
!
! - Geometry
!
    call jevech('PGEOMER', 'L', igeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)
!
! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
!
    call get_elas_id(zi(imate), elas_id)
!
! - Orthotropic parameters
!
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! - Compute RIGI_MECA
!
    do igau = 1, npg1
!
! ----- Compute matrix [B]: displacement -> strain (first order)
!
        call bmatmc(igau, nbsig, zr(igeom), ipoids, ivf, &
                    idfde, nno, nharm, jacgau, b)
!
! ----- Compute Hooke matrix [D]
!
        call dmatmc(fami, zi(imate), instan, '+', igau, &
                    1, angl_naut, nbsig, d)
!
! ----- Compute rigidity matrix [K] = [B]Tx[D]x[B]
!
        call btdbmc(b, d, jacgau, ndim, nno, &
                    nbsig, elas_id, btdb)
    end do
!
! - Set matrix in output field
!
    call jevech('PMATUUR', 'E', imatuu)
    k = 0
    do i = 1, nbinco
        do j = 1, i
            k = k+1
            zr(imatuu+k-1) = btdb(i, j)
        end do
    end do
!
end subroutine
