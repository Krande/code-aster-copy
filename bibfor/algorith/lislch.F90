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

subroutine lislch(lischa, ichar, charge)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnb.h"
    character(len=19) :: lischa
    integer(kind=8) :: ichar
    character(len=8) :: charge
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (LISTE_CHARGES)
!
! RETOURNE LE NOM DE LA CHARGE
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : SD LISTE DES CHARGES
! IN  ICHAR  : INDICE DE LA CHARGE
! OUT CHARGE : NOM DE LA CHARGE
!
!
!
!
    character(len=24) :: nomcha
    integer(kind=8) :: jncha
    integer(kind=8) :: nbchar
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    charge = ' '
    call lisnnb(lischa, nbchar)
!
    if (nbchar .ne. 0) then
        nomcha = lischa(1:19)//'.NCHA'
        call jeveuo(nomcha, 'L', jncha)
        charge = zk8(jncha-1+ichar)
    end if
!
    call jedema()
end subroutine
