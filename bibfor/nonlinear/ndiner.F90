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
subroutine ndiner(nbEqua, sddyna, hval_incr, hval_measse, cniner)
!
    use NonLinearDyna_module, only: compResiForce
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmdebg.h"
#include "asterfort/utmess.h"
#include "blas/dscal.h"
!
    integer(kind=8), intent(in) :: nbEqua
    character(len=19), intent(in) :: sddyna
    character(len=19), intent(in) :: hval_incr(*), hval_measse(*)
    character(len=19), intent(in) :: cniner
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute inertial force
!
! --------------------------------------------------------------------------------------------------
!
! In  nbEqua           : total number of equations
! In  sddyna           : datastructure for dynamic
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_measse      : hat-variable for matrix
! In  cniner           : name of field for inertia force
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: coefIner
    real(kind=8), pointer :: inerVale(:) => null()
    blas_int :: b_incx, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_32')
    end if
!
! - Coefficient
    coefIner = ndynre(sddyna, 'COEF_FORC_INER')
!
! - Compute
    call compResiForce(hval_incr, hval_measse, cniner)
    call jeveuo(cniner(1:19)//'.VALE', 'E', vr=inerVale)
    b_n = to_blas_int(nbEqua)
    b_incx = to_blas_int(1)
    call dscal(b_n, coefIner, inerVale, b_incx)
!
! - Debug
    if (niv .ge. 2) then
        call nmdebg('VECT', cniner, ifm)
    end if
!
end subroutine
