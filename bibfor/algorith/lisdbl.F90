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
subroutine lisdbl(lischa)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lislch.h"
#include "asterfort/lisnnb.h"
#include "asterfort/utmess.h"
    character(len=19) :: lischa
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (LISTE_CHARGES)
!
! VERIFICATIONS DES DOUBLONS
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : SD LISTE DES CHARGES
!
!
!
!
    integer(kind=8) :: ichar1, ichar2, nbchar
    character(len=8) :: charg1, charg2
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOMBRE DE CHARGES
!
    call lisnnb(lischa, nbchar)
!
! --- BOUCLE SUR LES CHARGES
!
    do ichar1 = 1, nbchar
        call lislch(lischa, ichar1, charg1)
        do ichar2 = 1, nbchar
            if (ichar1 .ne. ichar2) then
                call lislch(lischa, ichar2, charg2)
                if (charg1 .eq. charg2) then
                    call utmess('F', 'CHARGES5_2', sk=charg1)
                end if
            end if
        end do
    end do
!
    call jedema()
end subroutine
