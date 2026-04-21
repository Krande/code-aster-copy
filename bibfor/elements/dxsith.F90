! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine dxsith(nomte, materPara, sigma)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/dmatcp.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
!
    character(len=16), intent(in) :: nomte
    type(Material_Para), intent(inout) :: materPara
    real(kind=8), intent(out) :: sigma(*)
!
! --------------------------------------------------------------------------------------------------
!
!       CALCUL DES CONTRAINTES VRAIES
!        (==SIGMA_MECA - SIGMA_THER).
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbepsg = 8, nbcmp = 6
    character(len=8), parameter :: fami = "RIGI"
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: nbNode, npg
    integer(kind=8) :: i, j, icou, icpg, igauh, kpg, ipgh, iret, jvInstr, nbLayer
    integer(kind=8) :: npgh
    integer(kind=8) :: jnbspi, itab(8)
    real(kind=8) :: d(4, 4), time, epsth(nbepsg)
    aster_logical :: lDKTG
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, nno=nbNode, npg=npg)
!
    lDKTG = ASTER_FALSE
    if ((nomte .eq. 'MEDKTG3') .or. (nomte .eq. 'MEDKQG4')) then
        lDKTG = ASTER_TRUE
    end if

! - Get current time
    call tecach('ONO', 'PINSTR', 'L', iret, nval=8, itab=itab)
    jvInstr = itab(1)
    if (iret .eq. 0) then
        time = zr(jvInstr)
    else
        time = r8vide()
    end if

! - Get number of layers
    if (lDKTG) then
        nbLayer = 1
        npgh = 1
    else
        call jevech('PNBSP_I', 'L', jnbspi)
        npgh = 3
        nbLayer = zi(jnbspi-1+1)
        if (nbLayer .le. 0) then
            call utmess('F', 'ELEMENTS_46')
        end if
    end if

! - BOUCLE SUR LES POINTS DE GAUSS DE LA SURFACE
    do kpg = 1, npg
        do icou = 1, nbLayer
            do igauh = 1, npgh
                icpg = nbcmp*npgh*nbLayer*(kpg-1)+ &
                       nbcmp*npgh*(icou-1)+ &
                       nbcmp*(igauh-1)
                ipgh = npgh*(icou-1)+igauh

! ------------- Initializations of material parameters on current integration point
                call initParaPoin(kpg, igauh, materPara)

! ------------- Get thermal coefficient
                call verift('RIGI', kpg, ipgh, '+', materPara%jvMaterCode, epsth_=epsth(1))
                epsth(2) = epsth(1)
                epsth(3) = zero
                epsth(4) = zero
                epsth(5) = zero
                epsth(6) = zero

! ------------- Get Hooke matrix
                call dmatcp(materPara, "+", time, d)

! ------------- Compute stress
                do i = 1, 4
                    do j = 1, 4
                        sigma(icpg+i) = sigma(icpg+i)-epsth(j)*d(i, j)
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
