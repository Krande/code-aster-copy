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
subroutine ssrone(mag, isma, rota)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=8) :: mag, rota
    integer(kind=8) :: isma
! ----------------------------------------------------------------------
!     IN:  MAG : NOM DU MAILLAGE CONTENANT LA (SUPER)MAILLE ISMA
!          ISMA: NUMERO DE LA (SUPER)MAILLE DANS LE MAILLAGE MAG
!
!     OUT: ROTA: 'OUI' : LA ROTATION EST NECESSAIRE
!                'NON' : LA ROTATION N EST PAS NECESSAIRE
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    real(kind=8) :: r1
    integer(kind=8) :: rot1, rot2
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret, k
    real(kind=8), pointer :: para_r(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call jeexin(mag//'.PARA_R', iret)
    if (iret .gt. 0) then
!         -- ROT1= 1 : PEUT-ETRE , 0 : NON , 2 : OUI
        rot1 = 1
        call jeveuo(mag//'.PARA_R', 'L', vr=para_r)
    else
        rot1 = 0
    end if
    rot2 = rot1
    if (rot2 .eq. 1) then
        r1 = 0.0d0
        do k = 4, 6
            r1 = r1+abs(para_r(14*(isma-1)+k))
        end do
        rot1 = 0
        if (r1 .gt. 1.d-6) rot1 = 2
    end if
!
    ASSERT((rot1 .eq. 2) .or. (rot1 .eq. 0))
    if (rot1 .eq. 2) then
        rota = 'OUI'
    else if (rot1 .eq. 0) then
        rota = 'NON'
    end if
!
    call jedema()
end subroutine
