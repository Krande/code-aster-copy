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

subroutine ef0415(nomte)
!     CALCUL DE EFGE_ELNO
!     ------------------------------------------------------------------
    implicit none
#include "jeveux.h"
#include "asterfort/cosiro.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/utmess.h"
#include "asterfort/efcoq3d.h"

!
    character(len=16) :: nomte
!
!-----------------------------------------------------------------------
    integer(kind=8) ::  ichg
    integer(kind=8) ::  jcara, jeffg, jgeom
    integer(kind=8) :: lzi, lzr, nbcou
    integer(kind=8) :: npge, npgt
    integer(kind=8) :: nso

!-----------------------------------------------------------------------
    parameter(npge=3)
    parameter(npgt=10)
    integer(kind=8) ::  jmat, jnbspi
    integer(kind=8) :: nb1, nb2, npgsr, npgsn

!

    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    if (nomte .eq. 'MEC3QU9H') then
        nso = 4
    else if (nomte .eq. 'MEC3TR7H') then
        nso = 3
    end if
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PCACOQU', 'L', jcara)
!
!
    call cosiro(nomte, 'PCONTRR', 'L', 'UI', 'G', &
                ichg, 'S')
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
!
    if (nbcou .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if

    call jevete('&INEL.'//nomte//'.B', ' ', jmat)

    call jevech('PEFFORR', 'E', jeffg)

    call efcoq3d(nomte, nb1, nb2, zr(jcara), zr(jgeom), zr(lzr), &
                 zr(ichg), zr(jmat), zr(jeffg), &
                 nbcou, npgsn, npgsr, npge, nso, npgt)

!
end subroutine
