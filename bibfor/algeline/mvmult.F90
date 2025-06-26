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

subroutine mvmult(mat, vec, res)
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcmult.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/utmess.h"
#include "jeveux.h"

    character(len=*), intent(in) :: mat, vec, res

!    EFFECTUE LE PRODUIT D'UNE MATRICE PAR UN VECTEURS REELS. LE RESULTAT
!    EST STOCKE DANS N VECTEURS REELS
!     ATTENTION:
!       - MATRICE SYMETRIQUE OU NON, REELLE.
!       - LES VECTEURS INPUT ET OUTPUT REELS DOIVENT ETRE DISTINCTS
!       - POUR LES DDLS ELIMINES PAR AFFE_CHAR_CINE, ON NE PEUT PAS
!         CALCULER XSOL. CES DDLS SONT MIS A ZERO.
!     ---------------------------------------------------------------
!
    character(len=1) :: typmat, typres1, typres2
    character(len=19) :: vec19, res19, mat19
    integer(kind=8):: jchin, jchout, lmat
!
    call jemarq()

    mat19 = mat
    vec19 = vec
    res19 = res
!
    call mtdscr(mat19)
    call jeveuo(mat19//'.&INT', 'L', lmat)
    call jeveuo(vec19//'.VALE', 'L', jchin)
    call jeveuo(res19//'.VALE', 'E', jchout)
!
    call jelira(vec19//'.VALE', 'TYPE', cval=typres1)
    call jelira(res19//'.VALE', 'TYPE', cval=typres2)
!
    if (zi(lmat+3) .eq. 1) then
        typmat = 'R'
    else if (zi(lmat+3) .eq. 2) then
        typmat = 'C'
    else
        call utmess('F', 'ALGELINE2_86')
    end if
!
    ASSERT(typres1 == typres2)
    ASSERT(typres1 == typmat)
    !   PRODUIT MATRICE X VECTEUR :
!     ----------------------------------
    if (typres1 .eq. 'R') then
        call mrmult('ZERO', lmat, zr(jchin), zr(jchout), 1, ASTER_TRUE)
    else if (typres1 .eq. 'C') then
        call mcmult('ZERO', lmat, zc(jchin), zc(jchout), 1, ASTER_TRUE)
    end if

    call jedema()
end subroutine
