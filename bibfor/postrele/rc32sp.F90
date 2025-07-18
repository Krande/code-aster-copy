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
subroutine rc32sp(ze200, lieu, iocc1, iocc2, ns, &
                  sp, spmeca, instsp, nbsscyc, spss)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/rc32spa.h"
#include "asterfort/rc32spb.h"
!
    character(len=4) :: lieu
    integer(kind=8) :: iocc1, iocc2, ns, nbsscyc
    real(kind=8) :: sp(2), instsp(4), spmeca(2), spss(100)
    aster_logical :: ze200
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_ZE200
!     CALCUL DU SN
!
!     ------------------------------------------------------------------
! IN  : LIEU   : ='ORIG' : ORIGINE DU SEGEMNT, ='EXTR' : EXTREMITE
! OUT : SN     : PARTIE B3200 du SN
!
    integer(kind=8) :: n1, k
    character(len=8) :: methode
!
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
    do k = 1, 100
        spss(k) = 0.d0
    end do
    nbsscyc = 0
!
    call getvtx(' ', 'METHODE', scal=methode, nbret=n1)
    if (methode .eq. 'TRESCA') then
        call rc32spa(ze200, lieu, iocc1, iocc2, ns, &
                     sp, spmeca, instsp)
    else
        call rc32spb(ze200, lieu, iocc1, iocc2, ns, &
                     sp, spmeca, instsp, nbsscyc, spss)
    end if
!
    call jedema()
!
end subroutine
