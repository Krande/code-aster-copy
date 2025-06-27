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

subroutine nmcrpp(motfaz, iocc, prec, criter, tole)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
    character(len=*) :: motfaz
    integer(kind=8) :: iocc
    character(len=8) :: criter
    real(kind=8) :: prec, tole
!
! ----------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (UTILITAIRE - SELEC. INST.)
!
! LECTURE PRECISION/CRITERE
!
! ----------------------------------------------------------------------
!
! NB: SI LE CRITERE EST RELATIF MAIS QUE _PRECISION_ N'EST PAS
!     PRECISEE, ALORS PRECISION VAUT PREDEF
!
! IN  MOTFAC : MOT-FACTEUR POUR LIRE (LIST_INST/INST)
! IN  IOCC   : OCCURRENCE DU MOT-CLEF FACTEUR MOTFAC
! OUT PREC   : PRECISION DE RECHERCHE
! OUT CRITER : CRITERE DE SELECTION (RELATIF/ABSOLU)
! OUT TOLE   : TOLERANCE
!                +PREC POUR RELATIF
!                -PREC POUR ABSOLU

    integer(kind=8) :: n0, n1, n2
    character(len=16) :: motfac
    real(kind=8) :: predef
!
! ----------------------------------------------------------------------
!
    call jemarq()

    prec = 0.d0
    tole = 0.d0
    criter = 'RELATIF'
    motfac = motfaz
    predef = 1.d-6

!   CRITERE/PRECISION are only needed if INST or LIST_INST exist
    call getvr8(motfac, 'INST', iocc=iocc, nbret=n0)
    if (n0 .eq. 0) then
        call getvid(motfac, 'LIST_INST', iocc=iocc, nbret=n0)
    end if

    if (n0 .ne. 0) then
!       reading keywords
        call getvr8(motfac, 'PRECISION', iocc=iocc, scal=prec, nbret=n1)
        call getvtx(motfac, 'CRITERE', iocc=iocc, scal=criter, nbret=n2)
        if (criter .eq. 'ABSOLU') then
            if (n1 .eq. 0) then
                call utmess('F', 'LISTINST_1')
            end if
        else if (criter .eq. 'RELATIF') then
            if (n1 .eq. 0) then
                prec = predef
                call utmess('A', 'LISTINST_2', sr=predef)
            end if
        else
            ASSERT(.false.)
        end if

        if (prec .le. r8prem()) then
            call utmess('F', 'LISTINST_3')
        end if

!       setting 'tole'
        if (criter .eq. 'RELATIF') then
            tole = prec
        else if (criter .eq. 'ABSOLU') then
            tole = -prec
        else
            ASSERT(.false.)
        end if

    end if

    call jedema()

end subroutine
