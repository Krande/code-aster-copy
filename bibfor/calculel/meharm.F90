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

subroutine meharm(modele, nh, chharm)
    implicit none
#include "asterfort/dismoi.h"
#include "asterfort/mecact.h"
    character(len=*) :: modele
    character(len=24) :: chharm
    character(len=8) :: mailla
!
!    CETTE ROUTINE GENERE UN CHAMP D'HARMONIQUE (CARTE CONSTANTE)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nh
!-----------------------------------------------------------------------
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mailla)
    chharm = '&&MEHARM.NUME_HARM'
    call mecact('V', chharm, 'MAILLA', mailla, 'HARMON', &
                ncmp=1, nomcmp='NH', si=nh)
end subroutine
