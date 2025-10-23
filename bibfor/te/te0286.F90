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
subroutine te0286(option, nomte)
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ethdst.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/simtep.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in):: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D
!
! Options: EPOT_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: fami = "RIGI"
    real(kind=8) :: sigmEner(162), bsigmEner(81), anglNaut(3)
    real(kind=8) :: time, nharm
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: jvMater, jvGeom, jvDisp, jvHarm
    integer(kind=8) :: i, jvEner
    integer(kind=8) :: iret, ndim2
    integer(kind=8) :: nbsig, ndim, nno, npg
    real(kind=8) :: enerTherTher, enerPote
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
!
    ndim2 = 2
    if (lteatt('FOURIER', 'OUI')) then
        ndim = 3
    end if
    time = r8vide()
    nbsig = nbsigm()
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    call jevech('PMATERC', 'L', jvMater)

! - Orthotropic parameters
    call getElemOrientation(ndim2, nno, jvGeom, anglNaut)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get Fourier mode
    nharm = 0.d0
    call tecach('NNO', 'PHARMON', 'L', iret, iad=jvHarm)
    if (jvHarm .ne. 0) then
        nharm = dble(zi(jvHarm))
    end if

! - Compute "real" stress tensor at Gauss points
    call simtep(fami, nno, ndim, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), zr(jvDisp), &
                time, anglNaut, zi(jvMater), nharm, &
                sigmEner)

! - CALCUL DU VECTEUR DES FORCES INTERNES (BT*SIGMA)
    call bsigmc(nno, ndim, nbsig, npg, jvGaussWeight, &
                jvBaseFunc, jvDBaseFunc, zr(jvGeom), nharm, sigmEner, &
                bsigmEner)

! - CALCUL DU TERME EPSTH_T*D*EPSTH
    call ethdst(fami, nno, ndim, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), time, anglNaut, zi(jvMater), &
                enerTherTher)

! - CALCUL DE L'ENERGIE POTENTIELLE : 1/2*UT*K*U - UT*FTH + 1/2*EPSTHT*D*EPSTH
    enerPote = 0.d0
    do i = 1, ndim*nno
        enerPote = enerPote+bsigmEner(i)*zr(jvDisp+i-1)
    end do
    enerPote = enerPote+0.5d0*enerTherTher

! - Save energy
    call jevech('PENERDR', 'E', jvEner)
    zr(jvEner) = enerPote
!
end subroutine
