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

subroutine nmdecv(sddisc, nume_inst, i_event_acti, dtmin, retdec)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dinins.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
    character(len=19) :: sddisc
    integer(kind=8) :: nume_inst, i_event_acti, retdec
    real(kind=8) :: dtmin
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (GESTION DES EVENEMENTS - DECOUPE)
!
! VERIFICATIONS DE LA DECOUPE
!
! ----------------------------------------------------------------------
!
!
! In  sddisc           : datastructure for time discretization
! IN  NUMINS : NUMERO DE L'INSTANT COURANT
! In  i_event_acti     : index of active event
! IN  DTMIN  : INTERVALLE DE TEMPS MINIMAL SUR LA LISTE CREEE
! OUT RETDEC : CODE RETOUR DECOUPE
!               0 ECHEC DE LA DECOUPE
!               1 ON A DECOUPE
!               2 PAS DE DECOUPE
!
! ----------------------------------------------------------------------
!
    character(len=16) :: submet
    integer(kind=8) :: nbnivo, lenivo
    real(kind=8) :: pasmin
!
! ----------------------------------------------------------------------
!
!
!   Methode de decoupage
    call utdidt('L', sddisc, 'ECHE', 'SUBD_METHODE', index_=i_event_acti, &
                valk_=submet)
    ASSERT(submet .eq. 'MANUEL' .or. submet .eq. 'AUTO')

!
! --- PAS MINIMUM ATTEINT ?
!
    call utdidt('L', sddisc, 'ECHE', 'SUBD_PAS_MINI', index_=i_event_acti, &
                valr_=pasmin)

    if ((dtmin .lt. pasmin) .or. (dtmin .le. r8prem())) then
        retdec = 0
        call utmess('I', 'SUBDIVISE_16', sr=pasmin)
        goto 999
    else
        retdec = 1
    end if

!
! --- NIVEAU DE REDECOUPAGE MAXIMAL ATTEINT
!
    call utdidt('L', sddisc, 'ECHE', 'SUBD_NIVEAU', index_=i_event_acti, &
                vali_=nbnivo)

    ! Aucun controle par le nombre de decoupages dans la methode AUTO
    if (submet .eq. 'AUTO' .or. nbnivo .eq. -1) then
        retdec = 1
        goto 999
    end if

    ! Niveau courant (avant decoupe)
    lenivo = dinins(sddisc, nume_inst)

    if (lenivo-1 .ge. nbnivo) then
        call utmess('I', 'SUBDIVISE_17', si=nbnivo)
        retdec = 0
    else
        retdec = 1
    end if
!
999 continue

end subroutine
