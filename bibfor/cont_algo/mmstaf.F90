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

subroutine mmstaf(mesh, ndim, chdepd, coef_frot, &
                  nummae, aliase, nne, nummam, ksipc1, &
                  ksipc2, ksipr1, ksipr2, mult_lagr_f1, mult_lagr_f2, &
                  tang_1, tang_2, norm, pres_frot, dist_frot, &
                  indi_frot_eval)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mcopco.h"
#include "asterfort/mmvalp_scal.h"
#include "asterfort/mm_cycl_laugf.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: ndim
    character(len=19), intent(in) :: chdepd
    real(kind=8), intent(in) :: coef_frot
    integer(kind=8), intent(in) :: nummae
    character(len=8), intent(in) :: aliase
    integer(kind=8), intent(in) :: nne
    integer(kind=8), intent(in) :: nummam
    real(kind=8), intent(in) :: ksipc1
    real(kind=8), intent(in) :: ksipc2
    real(kind=8), intent(in) :: ksipr1
    real(kind=8), intent(in) :: ksipr2
    real(kind=8), intent(in) :: mult_lagr_f1(9)
    real(kind=8), intent(in) :: mult_lagr_f2(9)
    real(kind=8), intent(in) :: tang_1(3)
    real(kind=8), intent(in) :: tang_2(3)
    real(kind=8), intent(in) :: norm(3)
    real(kind=8), intent(out) :: pres_frot(3)
    real(kind=8), intent(out) :: dist_frot(3)
    integer(kind=8), intent(out) :: indi_frot_eval
!
! --------------------------------------------------------------------------------------------------
!
! Contact (continue method)
!
! Evaluate friction status
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ndim             : space size
! In  chdepd           : cumulated displacement
! In  coef_frot        : augmented ratio for friction
! In  lpenaf           : .true. if penalized friction
! In  nummae           : number of slave element
! In  aliase           : type of slave element
! In  nne              : number of nodes of slave element
! In  nummam           : number of master element
! In  kspic1           : first parametric coord. of contact point (in slave element)
! In  kspic2           : second parametric coord. of contact point (in slave element)
! In  kspir1           : first parametric coord. of projection of contact point (in master element)
! In  kspir2           : second parametric coord. of projection of contact point (in master element)
! In  mult_lagr_f1     : first lagrange multiplier for friction at nodes
! In  mult_lagr_f2     : second lagrange multiplier for friction at nodes
! In  tang_1           : first tangent vector
! In  tang_2           : second tangent vector
! In  norm             : normal
! Out pres_frot        : friction pressure
! Out dist_frot        : friction distance
! Out indi_frot_eval   : evaluation of new friction status
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ::  idim1, idim2
    real(kind=8) :: laug_frot_norm
    real(kind=8) :: dlagrf(2), dist_total(3)
    real(kind=8) :: ddeple(3), ddeplm(3)
    real(kind=8) :: mprojt(3, 3)
    real(kind=8) :: laug_frot(3)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    laug_frot_norm = 0.d0
    indi_frot_eval = 0
    laug_frot(:) = 0.d0
    dist_frot(:) = 0.d0
    dist_total(:) = 0.d0
    pres_frot(:) = 0.d0
    dlagrf(:) = 0.d0
    mprojt(:, :) = 0.d0
!
! - Tangent projection matrix
!
    do idim1 = 1, ndim
        do idim2 = 1, ndim
            mprojt(idim1, idim2) = -1.d0*norm(idim1)*norm(idim2)
        end do
    end do
    do idim1 = 1, ndim
        mprojt(idim1, idim1) = 1.d0+mprojt(idim1, idim1)
    end do
!
! - Lagrange multiplier for friction at current contact point
!
    call mmvalp_scal(ndim, aliase, nne, ksipc1, &
                     ksipc2, mult_lagr_f1, dlagrf(1))
    if (ndim .eq. 3) then
        call mmvalp_scal(ndim, aliase, nne, ksipc1, &
                         ksipc2, mult_lagr_f2, dlagrf(2))
    end if
!
! - Displacement increment
!
    call mcopco(mesh, chdepd, ndim, nummae, ksipc1, &
                ksipc2, ddeple)
    call mcopco(mesh, chdepd, ndim, nummam, ksipr1, &
                ksipr2, ddeplm)
!
! - Gap increment
!
    dist_total(1:3) = ddeple(1:3)-ddeplm(1:3)
!
! - Projection of gap increment on tangent plane
!
!    if (.not.lpenaf) then
    do idim1 = 1, ndim
        do idim2 = 1, ndim
            dist_frot(idim1) = mprojt(idim1, idim2)*dist_total(idim2)+dist_frot(idim1)
        end do
    end do
!    endif
!
! - Friction "pressure"
!
    if (ndim .eq. 2) then
        pres_frot(1:2) = dlagrf(1)*tang_1(1:2)
    else if (ndim .eq. 3) then
        pres_frot(1:3) = dlagrf(1)*tang_1(1:3)+dlagrf(2)*tang_2(1:3)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Norm of the augmented lagrangian for friction
!
    call mm_cycl_laugf(pres_frot, dist_frot, coef_frot, laug_frot_norm)
!
! - New status of friction (sign of augmented lagrangian)
!
    if (laug_frot_norm .le. 1.d0) then
        indi_frot_eval = 1
    else
        indi_frot_eval = 0
    end if
!
    call jedema()
end subroutine
