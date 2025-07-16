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

subroutine te0416(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/cosiro.h"
#include "asterfort/forngr.h"
#include "asterfort/fornpd.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
    character(len=16) :: option, nomte
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE : COQUE_3D
!     ----------------------------------------------------------------
!
!
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: ibid, iret, icompo
!
!
! DEB
!
    if (option .eq. 'FORC_NODA') then
!        -- PASSAGE DES CONTRAINTES DANS LE REPERE INTRINSEQUE :
        call cosiro(nomte, 'PSIEFR', 'L', 'UI', 'G', &
                    ibid, 'S')
    end if
!
!
    call tecach('ONO', 'PCOMPOR', 'L', iret, iad=icompo)
    if (icompo .eq. 0) then
        call fornpd(option, nomte)
        goto 999
    else

        call jevech('PCOMPOR', 'L', vk16=compor)

        if (compor(DEFO) .eq. 'GROT_GDEP') then

!           DEFORMATION DE GREEN

            call forngr(option, nomte)
!
            goto 999
!
        else if (compor(DEFO) (1:5) .eq. 'PETIT') then

            call fornpd(option, nomte)
        else

!----------- AUTRES MESURES DE DEFORMATIONS
!
            call utmess('F', 'ELEMENTS3_93', sk=compor(DEFO))
!
        end if
    end if
!
!
999 continue
!
!
!
end subroutine
