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

subroutine mepres(nomo, chpres, fonc, pres, cisa)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/mecact.h"
    character(len=8) :: nomo
    character(len=*) :: chpres
    real(kind=8) :: pres, cisa
!
    aster_logical :: fonc
!
!    CETTE ROUTINE GENERE UN CHAMP DE PRESSION CONSTANT (CARTE CONSTANTE)
!     FONC = .TRUE.  FORCE FONCTION
!     FONC = .FALSE. FORCE REELLE
!     PRES = COMPOSANTE DE PRESSION  DU CHAMP A GENERER
!     CISA = COMPOSANTE DE CISAILLEMENT  DU CHAMP A GENERER
!
!
    real(kind=8) :: rcmp(2)
!
    character(len=8) :: licmp(2), nomf(2), zero
    character(len=19) :: ligrmo
!
!
!-----------------------------------------------------------------------
    zero = '&FOZERO'
    ligrmo = nomo//'.MODELE    '
    licmp(1) = 'PRES'
    licmp(2) = 'CISA'
    rcmp(1) = pres
    rcmp(2) = cisa
    if (fonc) then
        nomf(1) = zero
        nomf(2) = zero
        call mecact('V', chpres, 'MODELE', ligrmo, 'PRES_F', &
                    ncmp=2, lnomcmp=licmp, vk=nomf)
    else
        call mecact('V', chpres, 'MODELE', ligrmo, 'PRES_R', &
                    ncmp=2, lnomcmp=licmp, vr=rcmp)
    end if
!
end subroutine
