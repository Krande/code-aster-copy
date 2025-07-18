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

function ndynin(sddyna, chaine)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
    integer(kind=8) :: ndynin
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=19) :: sddyna
    character(len=*) :: chaine
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (UTILITAIRE)
!
! INTERROGE SDDYNA POUR RENVOYER UN ENTIER
!
! ----------------------------------------------------------------------
!
!
! OUT NDYNIN : TYPE DE FORMULATION
!            -> 1 : FORMULATION EN DEPLACEMENT
!            -> 2 : FORMULATION EN VITESSE
!            -> 3 : FORMULATION EN ACCELERATION
! IN  SDDYNA : NOM DE LA SD DEDIEE A LA DYNAMIQUE
! IN  CHAINE :  = / 'FORMUL_CONTACT'
!                 / 'FORMUL_DYNAMIQUE'
!
!
!
!
!
    integer(kind=8) :: jtfor, jncha
    character(len=24) :: tfor, ncha
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ndynin = 0
!
! --- INTERROGATION
!
    tfor = sddyna(1:15)//'.TYPE_FOR'
    ncha = sddyna(1:15)//'.NBRE_CHA'
    call jeveuo(tfor, 'L', jtfor)
    call jeveuo(ncha, 'L', jncha)
    if (chaine .eq. 'FORMUL_DYNAMIQUE') then
        ndynin = zi(jtfor+1-1)
    else if (chaine .eq. 'NBRE_EXCIT') then
        ndynin = zi(jncha+1-1)
    else if (chaine .eq. 'NBRE_ONDE_PLANE') then
        ndynin = zi(jncha+2-1)
    else if (chaine .eq. 'NBRE_EXCIT_GENE') then
        ndynin = zi(jncha+3-1)
    else if (chaine .eq. 'NBRE_MODE_PROJ') then
        ndynin = zi(jncha+5-1)
    else
        ASSERT(.false.)
    end if
!
    call jedema()
end function
