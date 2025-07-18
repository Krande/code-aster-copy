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

subroutine cfnumm(defico, posnma, numnma)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=24), intent(in) :: defico
    integer(kind=8), intent(in) :: posnma
    integer(kind=8), intent(out) :: numnma
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES MAILLEES - UTILITAIRE)
!
! DONNE LES NUMEROS ABSOLUS DES MAILLES DE CONTACT
!
! ----------------------------------------------------------------------
!
!
! in  defico : sd de contact (definition)
! in  posnma : indice dans contno d'une maille
! out numnma : indice absolu de la maille dans le maillage
!
!
!
!
    character(len=24) :: contma
    integer(kind=8) :: jmaco
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- acces sd contact
!
    contma = defico(1:16)//'.MAILCO'
    call jeveuo(contma, 'L', jmaco)
!
! --- numero de la maille
!
    numnma = zi(jmaco+posnma-1)
!
    call jedema()
!
end subroutine
