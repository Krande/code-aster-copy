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
subroutine te0115(option, nomte)
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigvmc.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: AXIS_FOURIER
!
! Options: SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, npg, nbsig, i, dimmod
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: jvSigm, jvDisp, jvGeom, jvMater, jvHarm
    real(kind=8) :: sigm(54), anglNaut(3), time, nharm
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
    dimmod = 3
!
    nbsig = nbsigm()
    sigm = 0.d0
    ASSERT(nno .le. MT_NNOMAX2D)
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)

! - Current time
    time = r8vide()

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    call jevech('PMATERC', 'L', jvMater)

! - Orthotropic parameters:
    call getElemOrientation(ndim, nno, jvGeom, anglNaut)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get Fourier mode
    call jevech('PHARMON', 'L', jvHarm)
    nharm = dble(zi(jvHarm))

! - Compute mechanical stress (without effect of external state variables)
    call sigvmc('RIGI', nno, dimmod, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), zr(jvDisp), &
                time, anglNaut, zi(jvMater), nharm, sigm)

! - Final copy of stress
    call jevech('PCONTRR', 'E', jvSigm)
    do i = 1, nbsig*npg
        zr(jvSigm+i-1) = sigm(i)
    end do
!
end subroutine
