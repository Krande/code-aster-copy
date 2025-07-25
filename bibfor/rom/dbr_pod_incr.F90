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
subroutine dbr_pod_incr(lReuse, base, paraPod, q, s, &
                        v, nbModeOut, nbSnapOut)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dbr_calcpod_sele.h"
#include "asterfort/dbr_calcpod_svd.h"
#include "asterfort/detrsd.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/norm_frobenius.h"
#include "asterfort/romBaseCreate.h"
#include "asterfort/romTableRead.h"
#include "asterfort/romTableSave.h"
#include "asterfort/rsexch.h"
#include "asterfort/tbSuppressAllLines.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
#include "blas/dgesv.h"
!
    aster_logical, intent(in) :: lReuse
    type(ROM_DS_Empi), intent(in) :: base
    type(ROM_DS_ParaDBR_POD), intent(in) :: paraPod
    real(kind=8), pointer :: q(:), s(:), v(:)
    integer(kind=8), intent(out) :: nbModeOut, nbSnapOut
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE - Compute
!
! Incremental POD method
!
! --------------------------------------------------------------------------------------------------
!
! In  lReuse           : .true. if reuse
! In  base             : base
! In  paraPod          : datastructure for parameters (POD)
! Ptr q                : pointer to snapshots matrix (be modified after SVD)
! Ptr s                : pointer to singular values
! Ptr v                : pointer to singular vectors
! Out nbModeOut        : number of modes selected
! Out nbSnapOut        : number of snapshots used in incremental algorithm
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iAlgoIni, iAlgoEnd, iEqua, iSnap, iAlgoSnap, iAlgo, k, iMode, iCoorRedu
    integer(kind=8) :: nbEqua, nbSing, nbModeMaxi, nbSnapResult
    integer(kind=8) :: nbSnapPrev, nbModePrev
    real(kind=8) :: toleIncr, toleSVD
    character(len=8) :: baseName
    real(kind=8) :: norm_q, norm_r
    integer(kind=4) :: info
    real(kind=8), pointer :: qi(:) => null()
    real(kind=8), pointer :: ri(:) => null()
    real(kind=8), pointer :: rt(:) => null()
    real(kind=8), pointer :: vt(:) => null()
    real(kind=8), pointer :: g(:) => null()
    real(kind=8), pointer :: gt(:) => null()
    real(kind=8), pointer :: kv(:) => null()
    real(kind=8), pointer :: kt(:) => null()
    integer(kind=4), pointer :: IPIV(:) => null()
    real(kind=8), pointer :: b(:) => null()
    real(kind=8), pointer :: v_gamma(:) => null()
    character(len=24) :: mode, fieldName
    character(len=4) :: fieldSupp
    integer(kind=8) :: iret
    real(kind=8), pointer :: v_mode(:) => null()
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
    blas_int :: b_nrhs
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
!
! - Initializations
!
    nbModeOut = 0
    nbSnapOut = 0
!
! - Get parameters
!
    mode = '&&IPOD_MODE'
    nbEqua = base%mode%nbEqua
    nbModeMaxi = paraPod%nbModeMaxi
    nbSnapResult = paraPod%snap%nbSnap
    toleIncr = paraPod%toleIncr
    toleSVD = paraPod%toleSVD
    ASSERT(paraPod%baseType .eq. '3D')
!
! - Properties of previous base
!
    baseName = base%resultName
    nbModePrev = base%nbMode
    nbSnapPrev = base%nbSnap
    fieldName = base%mode%fieldName
    fieldSupp = base%mode%fieldSupp
!
! - Get previous reduced coordinates when reuse
!
    if (lReuse) then
        call romTableRead(paraPod%tablReduCoor)
        call tbSuppressAllLines(paraPod%tablReduCoor%tablResu%tablName)
    end if
!
! - Allocate objects
!
    AS_ALLOCATE(vr=qi, size=nbEqua)
    AS_ALLOCATE(vr=ri, size=nbEqua)
    AS_ALLOCATE(vr=rt, size=nbEqua)
    if (.not. lReuse) then
        ASSERT(nbModePrev .eq. 0)
        ASSERT(nbSnapPrev .eq. 0)
    end if
    AS_ALLOCATE(vr=vt, size=nbEqua*(nbSnapResult+nbModePrev))
    AS_ALLOCATE(vr=gt, size=(nbSnapResult+nbModePrev)*(nbSnapResult+nbSnapPrev))
    AS_ALLOCATE(vr=g, size=(nbSnapResult+nbModePrev)*(nbSnapResult+nbSnapPrev))
!
! - Initialize algorithm
!
    if (lReuse) then
! ----- Add previous modes in v
        do iMode = 1, nbModePrev
            call rsexch(' ', baseName, fieldName, iMode, mode, &
                        iret)
            if (fieldSupp .eq. 'NOEU') then
                call jeveuo(mode(1:19)//'.VALE', 'L', vr=v_mode)
            else if (fieldSupp .eq. 'ELGA') then
                call jeveuo(mode(1:19)//'.CELV', 'L', vr=v_mode)
            else
                ASSERT(ASTER_FALSE)
            end if
            do iEqua = 1, nbEqua
                vt(iEqua+nbEqua*(iMode-1)) = v_mode(iEqua)
            end do
        end do
! ----- Add previous reduced coordinates in gT
        do iCoorRedu = 1, nbModePrev*nbSnapPrev
            gt(iCoorRedu) = paraPod%tablReduCoor%coorRedu(iCoorRedu)
        end do
    else
! ----- Add first snap in v
        qi(1:nbEqua) = q(1:nbEqua)
        call norm_frobenius(nbEqua, qi, norm_q)
        if (norm_q .le. r8prem()) then
            norm_q = 1.d-16*sqrt(nbEqua*1.d0)
        end if
        vt(1:nbEqua) = qi(1:nbEqua)/norm_q
! ----- Add norm of first snap in gT
        gt(1) = norm_q
    end if
!
! - Suppress previous result datastructure
!
    if (lReuse) then
        call detrsd('RESULTAT', baseName)
        call romBaseCreate(base, nbModePrev)
    end if
!
! - Set bounds of algorithm
!
    if (lReuse) then
        iAlgoSnap = nbModePrev
        iAlgoIni = nbSnapPrev+1
        iAlgoEnd = nbSnapResult+nbSnapPrev
    else
        iAlgoSnap = 1
        iAlgoIni = 2
        iAlgoEnd = nbSnapResult
    end if
!
! - Main algorithm
!
    do iAlgo = iAlgoIni, iAlgoEnd
! ----- Get current snapshot
        if (lReuse) then
            do iEqua = 1, nbEqua
                qi(iEqua) = q(iEqua+nbEqua*(iAlgo-nbSnapPrev-1))
            end do
        else
            do iEqua = 1, nbEqua
                qi(iEqua) = q(iEqua+nbEqua*(iAlgo-1))
            end do
        end if
!
! ----- Compute norm of current snapshot
        call norm_frobenius(nbEqua, qi, norm_q)
        if (norm_q .le. r8prem()) then
            cycle
        end if
!
! ----- Compute {kt} = [v]^T {q} (projection of current snaphot on base)
        AS_ALLOCATE(vr=kt, size=iAlgoSnap)
        b_ldc = to_blas_int(iAlgoSnap)
        b_ldb = to_blas_int(nbEqua)
        b_lda = to_blas_int(nbEqua)
        b_m = to_blas_int(iAlgoSnap)
        b_n = to_blas_int(1)
        b_k = to_blas_int(nbEqua)
        call dgemm('T', 'N', b_m, b_n, b_k, &
                   1.d0, vt, b_lda, qi, b_ldb, &
                   0.d0, kt, b_ldc)
!
! ----- Compute [kv] = [v]^T [v]
        AS_ALLOCATE(vr=kv, size=iAlgoSnap*iAlgoSnap)
        b_ldc = to_blas_int(iAlgoSnap)
        b_ldb = to_blas_int(nbEqua)
        b_lda = to_blas_int(nbEqua)
        b_m = to_blas_int(iAlgoSnap)
        b_n = to_blas_int(iAlgoSnap)
        b_k = to_blas_int(nbEqua)
        call dgemm('T', 'N', b_m, b_n, b_k, &
                   1.d0, vt, b_lda, vt, b_ldb, &
                   0.d0, kv, b_ldc)
!
! ----- Solve [v]^T [v] {Y} = [v]^T {q} => {Y} are reduced coordinates
        AS_ALLOCATE(vi4=IPIV, size=iAlgoSnap)
        b_ldb = to_blas_int(iAlgoSnap)
        b_lda = to_blas_int(iAlgoSnap)
        b_n = to_blas_int(iAlgoSnap)
        b_nrhs = to_blas_int(1)
        call dgesv(b_n, b_nrhs, kv, b_lda, IPIV, &
                   kt, b_ldb, info)
!
! ----- Compute residu {r} = [v] {Y}
        b_ldc = to_blas_int(nbEqua)
        b_ldb = to_blas_int(iAlgoSnap)
        b_lda = to_blas_int(nbEqua)
        b_m = to_blas_int(nbEqua)
        b_n = to_blas_int(1)
        b_k = to_blas_int(iAlgoSnap)
        call dgemm('N', 'N', b_m, b_n, b_k, &
                   1.d0, vt, b_lda, kt, b_ldb, &
                   0.d0, rt, b_ldc)
        ri = qi-rt
!
! ----- Compute norm of residu
        call norm_frobenius(nbEqua, ri, norm_r)
!
! ----- Select vector or not ?
        if (norm_r/norm_q .ge. toleIncr) then
! --------- Add mode (residu !) at current iAlgoSnap iteration
            do iEqua = 1, nbEqua
                vt(iEqua+nbEqua*iAlgoSnap) = ri(iEqua)/norm_r
            end do
!
! --------- Add singular value
! --------- G matrice rectangulaire en quatre morceaux.
! --------- Partie du bas: iAlgoSnap+1
! --------- Partie du haut: 1 => iAlgoSnap
! --------- Partie gauche: 1 => iAlgo -1
! --------- Partie droit: iAlgo
!
! --------- Valeurs haut-gauche
            do iSnap = 1, iAlgoSnap
                g(iSnap+(iAlgoSnap+1)*(iAlgo-1)) = kt(iSnap)
                do k = 1, iAlgo-1
                    g(iSnap+(iAlgoSnap+1)*(k-1)) = gt(iSnap+iAlgoSnap*(k-1))
                end do
            end do
!
! --------- Valeurs bas-gauche
            do k = 1, iAlgo-1
                g((iAlgoSnap+1)*k) = 0.d0
            end do
!
! --------- Valeur bas-droite
            g((iAlgoSnap+1)*iAlgo) = norm_r
!
! --------- Valeurs haut-droite
            do k = 1, (iAlgoSnap+1)*iAlgo
                gt(k) = g(k)
            end do
!
! --------- Next snap
            iAlgoSnap = iAlgoSnap+1
!
        else
            do iSnap = 1, iAlgoSnap
                gt(iSnap+iAlgoSnap*(iAlgo-1)) = kt(iSnap)
            end do
        end if
        AS_DEALLOCATE(vr=kt)
        AS_DEALLOCATE(vr=kv)
        AS_DEALLOCATE(vi4=IPIV)
    end do
!
! - Deallocate objects
!
    AS_DEALLOCATE(vr=qi)
    AS_DEALLOCATE(vr=ri)
    AS_DEALLOCATE(vr=rt)
!
! - Final number of snapshots in base
!
    nbSnapOut = iAlgoEnd
!
! - Prepare matrix of reduced coordinates
!
    do iSnap = 1, iAlgoSnap*nbSnapOut
        g(iSnap) = gt(iSnap)
    end do
!
! - Compute SVD on matrix of reduced coordinates: Q = V S Wt
!
    call dbr_calcpod_svd(iAlgoSnap, nbSnapOut, g, s, b, &
                         nbSing)
!
! - Select empiric modes
!
    call dbr_calcpod_sele(nbModeMaxi, toleSVD, s, nbSing, nbModeOut)
!
! - Compute matrix of singular vector: V <= V * B (dim : [nbModeOut x nbEqua] )
!
    AS_ALLOCATE(vr=v, size=nbEqua*nbModeOut)
    b_ldc = to_blas_int(nbEqua)
    b_ldb = to_blas_int(iAlgoSnap)
    b_lda = to_blas_int(nbEqua)
    b_m = to_blas_int(nbEqua)
    b_n = to_blas_int(nbModeOut)
    b_k = to_blas_int(iAlgoSnap)
    call dgemm('N', 'N', b_m, b_n, b_k, &
               1.d0, vt, b_lda, b, b_ldb, &
               0.d0, v, b_ldc)
!
! - Compute reduced coordinates G <= B^T G (dim : [nbModeOut x nbSnapOut] )
!
    AS_ALLOCATE(vr=v_gamma, size=nbModeOut*nbSnapOut)
    b_ldc = to_blas_int(nbModeOut)
    b_ldb = to_blas_int(iAlgoSnap)
    b_lda = to_blas_int(iAlgoSnap)
    b_m = to_blas_int(nbModeOut)
    b_n = to_blas_int(nbSnapOut)
    b_k = to_blas_int(iAlgoSnap)
    call dgemm('T', 'N', b_m, b_n, b_k, &
               1.d0, b, b_lda, gt, b_ldb, &
               0.d0, v_gamma, b_ldc)
!
! - Save the reduced coordinates in a table
!
    if (niv .ge. 2) then
        call utmess('I', 'ROM5_39', ni=2, vali=[nbSnapOut, nbModeOut])
    end if
    do iSnap = 1, nbSnapOut
        call romTableSave(paraPod%tablReduCoor%tablResu, nbModeOut, v_gamma, numeSnap_=iSnap)
    end do
!
! - Debug print
!
    call utmess('I', 'ROM7_14', si=nbSnapOut)
!
! - Clean
!
    AS_DEALLOCATE(vr=v_gamma)
    AS_DEALLOCATE(vr=vt)
    AS_DEALLOCATE(vr=gt)
    AS_DEALLOCATE(vr=g)
    AS_DEALLOCATE(vr=b)
!
end subroutine
