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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine romAlgoNLCorrEFTherResidual(ds_algorom, vec2nd, cnvabt, cnresi, cn2mbr, &
                                       resi_rela, resi_maxi)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jeveuo.h"
!
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    character(len=24), intent(in) :: vec2nd, cnvabt, cnresi, cn2mbr
    real(kind=8), intent(out):: resi_rela, resi_maxi
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Solving non-linear problem THERMICS
!
! Evaluate residuals in applying HYPER-REDUCTION in CORR_EF phase
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_algorom       : datastructure for ROM parameters
! In  vec2nd           : applied loads
! In  cnvabt           : BT.T LAMBDA for Dirichlet loads
! In  cnresi           : non-linear residual
! In  cn2mbr           : equilibrium residual (to evaluate convergence)
! Out resi_rela        : value for RESI_GLOB_RELA
! Out resi_maxi        : value for RESI_GLOB_MAXI
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_hrom
    integer(kind=8) :: iEqua, nbEqua
    real(kind=8) :: vnorm
    real(kind=8), pointer :: v_cn2mbr(:) => null()
    real(kind=8), pointer :: v_cn2mbrr(:) => null()
    real(kind=8), pointer :: v_vec2nd(:) => null()
    real(kind=8), pointer :: v_cnvabt(:) => null()
    real(kind=8), pointer :: v_cnresi(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    resi_rela = 0.d0
    resi_maxi = 0.d0
    vnorm = 0.d0
!
! - Get parameters
!
    l_hrom = ds_algorom%l_hrom
    nbEqua = ds_algorom%ds_empi%mode%nbEqua
!
! - Access to vectors
!
    call jeveuo(cn2mbr(1:19)//'.VALE', 'E', vr=v_cn2mbr)
    call jeveuo(vec2nd(1:19)//'.VALE', 'L', vr=v_vec2nd)
    call jeveuo(cnvabt(1:19)//'.VALE', 'L', vr=v_cnvabt)
    call jeveuo(cnresi(1:19)//'.VALE', 'L', vr=v_cnresi)
!
! - Create residual
!
    do iEqua = 1, nbEqua
        v_cn2mbr(iEqua) = v_vec2nd(iEqua)-v_cnresi(iEqua)-v_cnvabt(iEqua)
    end do
!
! - Truncation of residual
!
    if (l_hrom) then
        do iEqua = 1, nbEqua
            if (ds_algorom%v_equa_sub(iEqua) .eq. 1) then
                v_vec2nd(iEqua) = 0.d0
                v_cnvabt(iEqua) = 0.d0
                v_cnresi(iEqua) = 0.d0
            end if
        end do
    end if
!
! - Product of modes by second member
!
    AS_ALLOCATE(vr=v_cn2mbrr, size=nbEqua)
!
! - Compute maximum
!
    do iEqua = 1, nbEqua
        v_cn2mbrr(iEqua) = v_vec2nd(iEqua)-v_cnresi(iEqua)-v_cnvabt(iEqua)
        resi_rela = resi_rela+(v_cn2mbrr(iEqua))**2
        vnorm = vnorm+(v_vec2nd(iEqua)-v_cnvabt(iEqua))**2
        resi_maxi = max(resi_maxi, abs(v_cn2mbrr(iEqua)))
    end do
!
! - Compute relative
!
    if (vnorm .gt. 0.d0) then
        resi_rela = sqrt(resi_rela/vnorm)
    end if
!
! - Cleaning
!
    AS_DEALLOCATE(vr=v_cn2mbrr)
!
end subroutine
