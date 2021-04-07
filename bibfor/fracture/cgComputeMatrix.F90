! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: tanguy.mathieu at edf.fr
!
subroutine cgComputeMatrix(cgField, cgTheta)
!
use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ismali.h"
#include "asterfort/gmatr1.h"
#include "asterfort/gmatc3.h"
#include "asterfort/imprsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"

!
    type(CalcG_field), intent(in)    :: cgField
    type(CalcG_theta), intent(inout) :: cgTheta
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Compute A Matrix from equation A*G(s)=g(theta) in 2D and 3D
!
!
! --------------------------------------------------------------------------------------------------
!
    integer          :: imatr, i, j
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    cgTheta%matrix='&&OP0060.MATRIX'

    if(cgField%ndim.eq.2)then
!       EN 2D, L'EQUATION EST SCALAIRE. A = 1     
        call wkvect(cgTheta%matrix, 'V V R8', 1, imatr)
        zr(imatr)=1.0
    else if(cgField%ndim.eq.3) then
!
        if(cgTheta%discretization .eq. 'LINEAIRE')then
!       Calcul de A dans le cas LINERAIRE
            call gmatc3(cgTheta%nb_fondNoeud,cgTheta%milieu,cgTheta%l_closed,&
                        cgTheta%absfond,cgTheta%matrix)
!
        elseif(cgTheta%discretization .eq. 'LEGENDRE') then
!!
            call wkvect(cgTheta%matrix, 'V V R8', (cgTheta%degree+1)*(cgTheta%degree+1), imatr)
            do i =1, cgTheta%degree+1
                do j=1, cgTheta%degree+1
                    if (i.eq.j) then
                        zr(imatr-1 +(i-1)*(cgTheta%degree+1)+j)=1.0
                    else
                        zr(imatr-1 +(i-1)*(cgTheta%degree+1)+j)=0.0
                    endif
                enddo
            enddo

        else
            ASSERT(.FALSE.)
        endif
    else
        ASSERT(.FALSE.)
    endif
    
    call jedema()
end subroutine
