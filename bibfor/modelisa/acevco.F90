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
subroutine acevco(nbocc, nlg, ier)
    implicit none
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbocc, nlg, ier
!                          AFFE_CARA_ELEM
!
!        VERIFICATION DES MOTS CLES POUR L'ELEMENT COQUE
!
! ----------------------------------------------------------------------
!  IN
!     NBOCC  : NOMBRE D'OCCURENCE
!  OUT
!     NLG    : NOMBRE TOTAL DE GROUPE DE MAILLE
!  IN/OUT
!     IER    : CUMUL DES ERREURS
! ----------------------------------------------------------------------
    integer(kind=8) :: ioc, nco, ne, nef, nex, nexf, ng, nin
    integer(kind=8) :: nk
    character(len=8) :: k8b, nomu
    character(len=16) :: concep, cmd
!-----------------------------------------------------------------------
    call getres(nomu, concep, cmd)
!
    nlg = 0
    do ioc = 1, nbocc
        call getvtx('COQUE', 'GROUP_MA', iocc=ioc, nbval=0, nbret=ng)
        call getvr8('COQUE', 'EPAIS', iocc=ioc, nbval=0, nbret=ne)
        call getvid('COQUE', 'EPAIS_FO', iocc=ioc, nbval=0, nbret=nef)
        call getvr8('COQUE', 'A_CIS', iocc=ioc, nbval=0, nbret=nk)
        call getvr8('COQUE', 'EXCENTREMENT', iocc=ioc, nbval=0, nbret=nex)
        call getvid('COQUE', 'EXCENTREMENT_FO', iocc=ioc, nbval=0, nbret=nexf)
        call getvtx('COQUE', 'INER_ROTA', iocc=ioc, nbval=0, nbret=nin)
        call getvtx('COQUE', 'MODI_METRIQUE', iocc=ioc, nbval=0, nbret=nco)
!
        if (ioc .eq. 1 .and. abs(ne+nef) .ne. 1) then
            call utmess('E', 'MODELISA_53')
            ier = ier+1
        end if
!
        if ((nex+nexf) .ne. 0 .and. nin .ne. 0) then
            call getvtx('COQUE', 'INER_ROTA', iocc=ioc, scal=k8b, nbret=nin)
            if (k8b .eq. 'NON') then
                call utmess('E', 'MODELISA_54')
                ier = ier+1
            end if
        end if
!
        nlg = max(nlg, -ng)
    end do
!
end subroutine
