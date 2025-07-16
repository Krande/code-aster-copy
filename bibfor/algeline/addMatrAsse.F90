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

subroutine addMatrAsse(mat1, mat2, coeff1, coeff2, matres)
    implicit none
!
!     COMBINAISON LINEAIRE DE MATRICE
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/utmess.h"
#include "asterfort/vrrefe.h"
!
    character(len=*), intent(in) :: mat1, mat2, matres
    real(kind=8), intent(in) :: coeff1, coeff2

!
    integer(kind=8), parameter :: nbmatr = 2
    character(len=19), parameter :: nameMatr = "&MATRADD"
    character(len=19) :: listMatr(2), mat19
    integer(kind=8) :: iocc, jrefe, jpomr, iret
    aster_logical :: reuse
! ------------------------------------------------------------------
!
    call jemarq()
!
    reuse = ASTER_FALSE
    if (matres == mat1) then
        mat19 = nameMatr
        reuse = ASTER_TRUE
    elseif (matres == mat2) then
        mat19 = nameMatr
        reuse = ASTER_TRUE
    else
        mat19 = matres
    end if

    listMatr(1:2) = [mat1, mat2]
    jpomr = 1
    do iocc = 1, nbMatr
!       -- on recherche une eventuelle matrice non symetrique
        call jeveuo(listMatr(iocc)//'.REFA', 'L', jrefe)
        if (zk24(jrefe-1+9) .eq. 'MR') then
            jpomr = iocc
        end if
    end do
!
!   --- controle des references :
!   --------------------------------
    call vrrefe(listMatr(1), listMatr(2), iret)
    if (iret .ne. 0) then
        call utmess('F', 'ALGELINE2_28', nk=2, valk=listMatr)
    end if
!
!   -- combinaison des matrices :
!   ------------------------------------------------------------------
! initialisation de la matrice resultat :
    call mtdefs(mat19, listMatr(jpomr), 'G', ' ')
    call mtcmbl(nbMatr, ['R', 'R'], [coeff1, coeff2], listMatr, mat19, &
                ' ', ' ', 'ELIM=')
!
    if (reuse) then
        call detrsd('MATR_ASSE', matres)
        call copisd('MATR_ASSE', 'G', mat19, matres)
        call detrsd('MATR_ASSE', mat19)
    end if

    call jedema()
end subroutine
