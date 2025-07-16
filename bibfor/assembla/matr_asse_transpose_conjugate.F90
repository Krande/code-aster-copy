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
!
subroutine matr_asse_transpose_conjugate(matas)
! person_in_charge: nicolas.tardieu at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
!
    character(len=*) :: matas
!-----------------------------------------------------------------------
! But : transposer une matrice assemblee
!
! in/jxvar  k* matas : nom de la sd_matr_asse a transposer
!
!-----------------------------------------------------------------------
    character(len=19) :: matas1
    character(len=3) :: tysca
    integer(kind=8) :: n1, neq, jvalm1, jvalm2, jvaltmp, i
    blas_int :: b_incx, b_incy, b_n
!-------------------------------------------------------------------
    call jemarq()
    matas1 = matas
! symetrie
    call jelira(matas1//'.VALM', 'NUTIOC', n1)
    ASSERT(n1 .eq. 1 .or. n1 .eq. 2)
! reelle ou complexe
    call jelira(matas1//'.VALM', 'TYPE', cval=tysca)
! taille des vecteurs
    call jelira(jexnum(matas1//'.VALM', 1), 'LONMAX', neq)
! matrice symetrique
! ------------------
    if (n1 .eq. 1) then
        if (tysca .eq. 'R') then
            goto 999
        else if (tysca .eq. 'C') then
            call jeveuo(jexnum(matas1//'.VALM', 1), 'E', jvalm1)
            do i = 1, neq
                zc(jvalm1-1+i) = dconjg(zc(jvalm1-1+i))
            end do
        else
            ASSERT(.false.)
        end if
    else
! matrice non-symetrique
! ----------------------
        call wkvect('&&matr_transpose.VALM', 'V V '//tysca, neq, jvaltmp)
        call jeveuo(jexnum(matas1//'.VALM', 1), 'E', jvalm1)
        call jeveuo(jexnum(matas1//'.VALM', 2), 'E', jvalm2)
        if (tysca .eq. 'R') then
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvalm1), b_incx, zr(jvaltmp), b_incy)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvalm2), b_incx, zr(jvalm1), b_incy)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvaltmp), b_incx, zr(jvalm2), b_incy)
        else if (tysca .eq. 'C') then
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, zc(jvalm1), b_incx, zc(jvaltmp), b_incy)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, zc(jvalm2), b_incx, zc(jvalm1), b_incy)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, zc(jvaltmp), b_incx, zc(jvalm2), b_incy)
            do i = 1, neq
                zc(jvalm1-1+i) = dconjg(zc(jvalm1-1+i))
                zc(jvalm2-1+i) = dconjg(zc(jvalm2-1+i))
            end do
        else
            ASSERT(.false.)
        end if
        call jedetr('&&matr_transpose.VALM')
    end if
!
999 continue
    call jedema()
end subroutine
