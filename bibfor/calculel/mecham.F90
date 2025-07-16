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
subroutine mecham(optionZ, modelZ, cara, nh, chgeoz, &
                  chcara, chharz, iret)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/utmess.h"
#include "asterfort/checkSuperElement.h"
!
    integer(kind=8) :: iret, nh
    character(len=*), intent(in) :: optionZ, modelZ
    character(len=*) :: cara
    character(len=*) :: chgeoz, chcara(*), chharz
!
! --------------------------------------------------------------------------------------------------
!
!     VERIFICATION DES CHAMPS DONNES :
!        - ON VERIFIE LES EVENTUELLES SOUS-STRUCTURES STATIQUES
!        - ON VERIFIE S'IL Y A 1 LIGREL DANS LE MODELE
!
! --------------------------------------------------------------------------------------------------
!
! IN  : OPTION : OPTION DE CALCUL
! IN  : MODELE : MODELE
! IN  : CARA   : CHAMP DE CARA_ELEM
! IN  : NH     : NUMERO D'HARMONIQUE DE FOURIER
! OUT : CHGEOZ : NOM DE CHAMP DE GEOMETRIE TROUVE
! OUT : CHCARA : NOMS DES CHAMPS DE CARACTERISTIQUES TROUVES
! OUT : CHHARZ : NOM DU CHAMP D'HARMONIQUE DE FOURIER TROUVE
! OUT : IRET  : CODE RETOUR
!                = 0 : LE MODELE CONTIENT DES ELEMENTS FINIS
!                = 1 : LE MODELE NE CONTIENT PAS D'ELEMENTS FINIS
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model, exiele
    character(len=24) :: chgeom, chharm
    character(len=16) :: option
    integer(kind=8) :: nbSuperElement
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    chgeom = ' '
    chharm = ' '
    option = optionZ
!
    ASSERT(modelZ(1:1) .ne. ' ')
    model = modelZ

! - Check if super-elements have been computed
    call checkSuperElement(option, model)
!
!
!     --- ON REGARDE S'IL Y A 1 LIGREL DANS LE MODELE ---
    call dismoi('EXI_ELEM', model, 'MODELE', repk=exiele)
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSuperElement)
    if (exiele(1:3) .eq. 'OUI') then
        iret = 0
    else
        iret = 1
    end if
    if (iret .eq. 1 .and. nbSuperElement .eq. 0) then
        call utmess('F', 'CALCULEL3_35')
    end if
!
    if (iret .ne. 1) then
        call megeom(model, chgeom)
        call mecara(cara, chcara)
        call meharm(model, nh, chharm)
    end if
!
    chgeoz = chgeom
    chharz = chharm
    call jedema()
end subroutine
