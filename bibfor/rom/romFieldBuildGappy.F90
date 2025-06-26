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
subroutine romFieldBuildGappy(resultRom, fieldBuild)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
#include "blas/dgesv.h"
!
    type(ROM_DS_Result), intent(in) :: resultRom
    type(ROM_DS_FieldBuild), intent(inout) :: fieldBuild
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Field build
!
! Compute reduced coordinates with Gappy-POD
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
    type(ROM_DS_Field) :: mode
    character(len=8) :: resultRomName
    integer(kind=8) :: nbStore, nbEqua, nbEquaRom, nbMode, nbEquaRID
    integer(kind=8) :: iStore, iMode, iEqua
    integer(kind=8) :: numeStore, numeEqua, iret
    integer(kind=4) :: systInfo
    real(kind=8), pointer :: systMatr(:) => null(), systVect(:) => null()
    integer(kind=4), pointer :: systPerm(:) => null()
    character(len=4) :: fieldSupp
    character(len=24) :: fieldRom, fieldName
    real(kind=8), pointer :: valeField(:) => null(), valeRom(:) => null()
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
    blas_int :: b_nrhs
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM17_8')
    end if
!
! - Get parameters about reduced results
!
    resultRomName = resultRom%resultName
    nbStore = resultRom%nbStore
!
! - Get parameters about base
!
    base = fieldBuild%base
    mode = base%mode
    nbMode = base%nbMode
!
! - Get parameters about mode
!
    fieldName = mode%fieldName
!
! - Get parameters
!
    nbEquaRID = fieldBuild%nbEquaRID
!
! - Allocate objects
!
    AS_ALLOCATE(vr=fieldBuild%reduMatr, size=nbMode*(nbStore-1))
    AS_ALLOCATE(vr=systMatr, size=nbMode*nbMode)
    AS_ALLOCATE(vr=systVect, size=nbMode)
    AS_ALLOCATE(vi4=systPerm, size=nbMode)
    if (fieldBuild%lRIDTrunc) then
        AS_ALLOCATE(vr=valeField, size=nbEquaRID)
    end if
!
! - Compute Gappy POD
!
    do iStore = 1, nbStore-1
        numeStore = iStore
!
! ----- Get current field from reduced model
        call rsexch(' ', resultRomName, fieldName, numeStore, fieldRom, &
                    iret)
        ASSERT(iret .eq. 0)
        call dismoi('TYPE_CHAMP', fieldRom, 'CHAMP', repk=fieldSupp)
        if (fieldSupp == 'NOEU') then
            call jeveuo(fieldRom(1:19)//'.VALE', 'L', vr=valeRom)
            call jelira(fieldRom(1:19)//'.VALE', 'LONMAX', nbEquaRom)
        else if (fieldSupp == 'ELGA') then
            call jeveuo(fieldRom(1:19)//'.CELV', 'L', vr=valeRom)
            call jelira(fieldRom(1:19)//'.CELV', 'LONMAX', nbEquaRom)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ----- Truncate input field if required
        if (fieldBuild%lRIDTrunc) then
            nbEqua = 0
            do iEqua = 1, nbEquaRom
                numeEqua = fieldBuild%equaRIDTrunc(iEqua)
                if (numeEqua .ne. 0) then
                    valeField(numeEqua) = valeRom(iEqua)
                    nbEqua = nbEqua+1
                end if
            end do
            ASSERT(nbEqua .eq. nbEquaRID)
        else
            ASSERT(nbEquaRom .eq. nbEquaRID)
            if (fieldSupp == 'NOEU') then
                call jeveuo(fieldRom(1:19)//'.VALE', 'L', vr=valeField)
            else if (fieldSupp == 'ELGA') then
                call jeveuo(fieldRom(1:19)//'.CELV', 'L', vr=valeField)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
! ----- Compute matrix and vector
        b_ldc = to_blas_int(nbMode)
        b_ldb = to_blas_int(nbEquaRID)
        b_lda = to_blas_int(nbEquaRID)
        b_m = to_blas_int(nbMode)
        b_n = to_blas_int(1)
        b_k = to_blas_int(nbEquaRID)
        call dgemm('T', 'N', b_m, b_n, b_k, &
                   1.d0, fieldBuild%matrPhiRID, b_lda, valeField, b_ldb, &
                   0.d0, systVect, b_ldc)
        b_ldc = to_blas_int(nbMode)
        b_ldb = to_blas_int(nbEquaRID)
        b_lda = to_blas_int(nbEquaRID)
        b_m = to_blas_int(nbMode)
        b_n = to_blas_int(nbMode)
        b_k = to_blas_int(nbEquaRID)
        call dgemm('T', 'N', b_m, b_n, b_k, &
                   1.d0, fieldBuild%matrPhiRID, b_lda, fieldBuild%matrPhiRID, b_ldb, &
                   0.d0, systMatr, b_ldc)
!
! ----- Solve system
        b_ldb = to_blas_int(nbMode)
        b_lda = to_blas_int(nbMode)
        b_n = to_blas_int(nbMode)
        b_nrhs = to_blas_int(1)
        call dgesv(b_n, b_nrhs, systMatr, b_lda, systPerm, &
                   systVect, b_ldb, systInfo)
        if (systInfo .ne. 0) then
            call utmess('F', 'ROM17_9')
        end if
!
! ----- Copy result
        do iMode = 1, nbMode
            fieldBuild%reduMatr(iMode+nbMode*(numeStore-1)) = systVect(iMode)
        end do
    end do
!
! - Clean
!
    AS_DEALLOCATE(vr=systMatr)
    AS_DEALLOCATE(vr=systVect)
    AS_DEALLOCATE(vi4=systPerm)
    if (fieldBuild%lRIDTrunc) then
        AS_DEALLOCATE(vr=valeField)
    end if
!
end subroutine
