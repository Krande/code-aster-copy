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

subroutine cfnomm(noma, defico, typent, posent, noment)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfnumm.h"
#include "asterfort/cfnumn.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: noma
    character(len=24), intent(in) :: defico
    integer(kind=8), intent(in) :: posent
    character(len=4), intent(in) :: typent
    character(len=8), intent(out) :: noment
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES MAILLEES - UTILITAIRE)
!
! DONNE LE NOM DE L'ENTITE A PARTIR DE SON NUMERO
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  DEFICO : SD DE DEFINITION DU CONTACT (ISSUE D'AFFE_CHAR_MECA)
! IN  POSENT : POSITION DE L'ENTITE DANS LES SD CONTACT
! IN  TYPENT : TYPE D'ENTITE
!                <MAIL>  MAILLE
!                <NOEU>  NOEUD
! OUT NOMENT : NOM DE L'ENTITE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nummai, numnoe(1)
    integer(kind=8) :: posmai, posnoe(1)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    noment = ' '
!
! --- PREPARATION DES CHAINES POUR LES NOMS
!
    if (typent .eq. 'MAIL') then
        posmai = posent
        call cfnumm(defico, posmai, nummai)
        noment = int_to_char8(nummai)
!
    else if (typent .eq. 'NOEU') then
        posnoe = abs(posent)
        call cfnumn(defico, 1, posnoe(1), numnoe(1))
        noment = int_to_char8(numnoe(1))
    else
        ASSERT(.false.)
    end if
!
    call jedema()
!
end subroutine
