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
subroutine cuelpv(numlia, resocu, nbliai, lelpiv)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jeveuo.h"
    aster_logical :: lelpiv
    integer(kind=8) :: numlia
    integer(kind=8) :: nbliai
    character(len=24) :: resocu
!
! ----------------------------------------------------------------------
! ROUTINE APPELEE PAR : ALGOCU
! ----------------------------------------------------------------------
!
!  SAVOIR SI LA LIAISON EST A PIVOT NUL
!
! IN  RESOCU : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NUMLIA : NUMERO DE LA LIAISON
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
! OUT CFELPV : .TRUE.  SI LIAISON A PIVOT NUL
!              .FALSE. SINON
!
!
!
!
    integer(kind=8) :: iote
    character(len=19) :: liot
    integer(kind=8) :: jliot
! ======================================================================
    liot = resocu(1:14)//'.LIOT'
    call jeveuo(liot, 'L', jliot)
! ======================================================================
    lelpiv = .false.
    do iote = 1, zi(jliot+4*nbliai)
        if (zi(jliot-1+iote) .eq. numlia) then
            lelpiv = .true.
            goto 100
        end if
    end do
!
100 continue
! ======================================================================
end subroutine
