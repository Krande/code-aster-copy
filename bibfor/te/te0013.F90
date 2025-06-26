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

subroutine te0013(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/metau1.h"
#include "asterfort/metau2.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/sigtmc.h"
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
! Elements: 2D et 3D
! Option: CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4) :: fami
    real(kind=8) :: bsigma(81), sigth(162), angl_naut(3), time, nharm
    integer(kind=8) :: i, idfde, igeom, imate, ipoids, itemps, ivectu, iret
    integer(kind=8) :: ivf, nbsig, ndim, nno, npg
    real(kind=8) :: zero
    aster_logical :: l_meta
!
! --------------------------------------------------------------------------------------------------
!
    zero = 0.d0
    time = r8vide()
    nharm = zero
    ndim = 2
    fami = 'RIGI'
    sigth(:) = zero
    bsigma(:) = zero
!
!
! - Finite element informations
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, jpoids=ipoids, &
                     jvf=ivf, jdfde=idfde)
!

! - Compute CHAR_MECA_TEMP_R for metallurgy
!
    if (ndim .eq. 3) then
        call metau2(l_meta)
    else
        call metau1(l_meta)
    end if

    if (l_meta) then
        goto 40
    end if

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
! - Orthotropic parameters
!
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! - Get time
!
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) then
        time = zr(itemps)
    end if
!
! - Compute thermal stresses {SIGTH}
!
    call sigtmc('RIGI', ndim, nbsig, npg, &
                time, zi(imate), angl_naut, &
                option, sigth)
!
! - Compute CHAR_MECA_TEMP_R: [B]Tx{SIGTH}
!
    call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                ivf, idfde, zr(igeom), nharm, sigth, &
                bsigma)
!
! - Set output vector
!
    call jevech('PVECTUR', 'E', ivectu)
!
    do i = 1, ndim*nno
        zr(ivectu+i-1) = bsigma(i)
    end do
!
40  continue
end subroutine
