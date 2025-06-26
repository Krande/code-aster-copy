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
! person_in_charge: sylvie.granet at edf.fr
!
subroutine calc2nd(ds_thm, j_mater, &
                   lMatr, lSigm, &
                   ndim, dimdef, dimcon, &
                   adde2nd, adco2nd, &
                   defgem, defgep, &
                   congem, congep, &
                   dsde)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dil2gr.h"
#include "asterfort/dilpen.h"

    type(THM_DS), intent(inout) :: ds_thm
    aster_logical, intent(in) :: lMatr, lSigm
    integer(kind=8), intent(in) :: j_mater, ndim, dimdef, dimcon
    integer(kind=8), intent(in) :: adde2nd, adco2nd
    real(kind=8), intent(in) :: defgem(dimdef)
    real(kind=8), intent(in) :: defgep(dimdef)
    real(kind=8), intent(in) :: congem(dimcon)
    real(kind=8), intent(inout) :: congep(dimcon)
    real(kind=8), intent(inout) :: dsde(dimcon, dimdef)
!
! --------------------------------------------------------------------------------------------------
!
! THM - Second gradient
!
! Compute generalized stresses and matrix for second gradient
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  j_mater          : coded material address
! In  ndim             : dimension of space (2 or 3)
! In  dimdef           : dimension of generalized strains vector
! In  dimcon           : dimension of generalized stresses vector
! In  adde2nd          : adress of second gradient in generalized strains vector
! In  adco2nd          : adress of second gradient in generalized stresses vector
! In  defgem           : generalized strains - At begin of current step
! In  defgep           : generalized strains - At end of current step
! In  congem           : generalized stresses - At begin of current step
! IO  congep           : generalized stresses - At end of current step
! Out dsde             : derivative matrix stress/strain (behaviour only)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_dim
    real(kind=8) :: rpena
    real(kind=8) :: eplcp(ndim), silcp(ndim)
    real(kind=8) :: dsde2g(ndim, ndim)

! - Get penalization coefficient
    call dilpen(j_mater, rpena)

! - Get second gradient behaviour generalized strains at end of current steps
    eplcp(1:ndim) = defgep(adde2nd+2:adde2nd+1+ndim)
    silcp(:) = 0.d0
    dsde2g(:, :) = 0.d0

! - Compute second gradient behaviour
    call dil2gr(j_mater, ndim, ndim, eplcp, silcp, dsde2g)

! - Compute second gradient generalized stresses at end of current step
    if (lSigm) then
        congep(adco2nd) = defgep(adde2nd+2+ndim)+ &
                          rpena*(defgep(adde2nd)-defgep(adde2nd+1))
        congep(adco2nd+1) = -defgep(adde2nd+2+ndim)- &
                            rpena*(defgep(adde2nd)-defgep(adde2nd+1))
        do i_dim = 1, ndim
            congep(adco2nd+1+i_dim) = silcp(i_dim)
        end do
        congep(adco2nd+2+ndim) = defgep(adde2nd)-defgep(adde2nd+1)
    end if

! - Compute second gradient matrix
    if (lMatr) then
        dsde(adco2nd, adde2nd) = rpena
        dsde(adco2nd, adde2nd+1) = -rpena
        dsde(adco2nd, adde2nd+2+ndim) = 1.0d0
        dsde(adco2nd+1, adde2nd) = -rpena
        dsde(adco2nd+1, adde2nd+1) = rpena
        dsde(adco2nd+1, adde2nd+2+ndim) = -1.0d0
        dsde(adco2nd+2:adco2nd+1+ndim, adde2nd+2:adde2nd+1+ndim) = dsde2g(1:ndim, 1:ndim)
        dsde(adco2nd+2+ndim, adde2nd) = 1.0d0
        dsde(adco2nd+2+ndim, adde2nd+1) = -1.0d0
    end if

end subroutine
