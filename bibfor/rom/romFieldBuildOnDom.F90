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
subroutine romFieldBuildOnDom(resultRom, fieldBuild)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
!
    type(ROM_DS_Result), intent(in) :: resultRom
    type(ROM_DS_FieldBuild), intent(inout) :: fieldBuild
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Field build
!
! Construct field on complete domain and all storing index
!
! --------------------------------------------------------------------------------------------------
!
! In  base             : base
! In  resultRom        : reduced results
! IO  fieldBuild       : field to reconstruct
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    type(ROM_DS_Empi) :: base
    integer(kind=8) :: nbMode, nbEqua, nbStore
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM17_6')
    end if
!
! - Get parameters
!
    base = fieldBuild%base
    nbMode = base%nbMode
    nbEqua = base%mode%nbEqua
    nbStore = resultRom%nbStore
!
! - Allocate object
!
    AS_ALLOCATE(vr=fieldBuild%fieldTransientVale, size=nbEqua*(nbStore-1))
!
! - Construct object
!
    b_ldc = to_blas_int(nbEqua)
    b_ldb = to_blas_int(nbMode)
    b_lda = to_blas_int(nbEqua)
    b_m = to_blas_int(nbEqua)
    b_n = to_blas_int(nbStore-1)
    b_k = to_blas_int(nbMode)
    call dgemm('N', 'N', b_m, b_n, b_k, &
               1.d0, fieldBuild%matrPhi, b_lda, fieldBuild%reduMatr, b_ldb, &
               0.d0, fieldBuild%fieldTransientVale, b_ldc)
!
end subroutine
