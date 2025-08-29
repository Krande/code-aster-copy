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

subroutine cclpco(option, &
                  resultOut, numeStore, &
                  nbParaOut, lpaout, lchout)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: resultOut
    integer(kind=8), intent(in) :: numeStore
    integer(kind=8), intent(out) :: nbParaOut
    character(len=8), intent(out) :: lpaout(1)
    character(len=24), intent(out) :: lchout(1)
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! Compute ELEM, ELNO and ELGA fields - Set output fields
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION A CALCULER
!   RESUOU  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT OUT
!   NUMORD  I    NUMERO D'ORDRE COURANT
!
! OUT :
!   NBPAOU  I    NOMBRE DE PARAMETRES OUT
!   LIPAOU  K8*  LISTE DES PARAMETRES OUT
!   LICHOU  K8*  LISTE DES CHAMPS OUT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: optionNume, cataNbIn, ipara, ierd
    integer(kind=8) :: cataNbOut, kpara, physQuanNume
    character(len=8) :: physQuanName, scalType
    character(len=19) :: fieldOut
    integer(kind=8), pointer :: cataDescopt(:) => null()
    character(len=8), pointer :: cataOptpara(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nbParaOut = 0
    lpaout = " "
    lchout = " "

! - Access to catalog of options
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), optionNume)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', optionNume), 'L', vi=cataDescopt)
    call jeveuo(jexnum('&CATA.OP.OPTPARA', optionNume), 'L', vk8=cataOptpara)
    cataNbIn = cataDescopt(2)
    cataNbOut = cataDescopt(3)

! - Look for real
    if (cataNbOut .eq. 1) then
        ipara = 1

    elseif (cataNbOut .eq. 2) then
        ipara = 0
        do kpara = 1, 2
            physQuanNume = cataDescopt(4+cataNbIn+kpara)
            call jenuno(jexnum('&CATA.GD.NOMGD', physQuanNume), physQuanName)
            call dismoi('TYPE_SCA', physQuanName, 'GRANDEUR', repk=scalType)
            if (scalType .eq. 'R') ipara = kpara
        end do
        ASSERT(ipara .gt. 0)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Set output field
    nbParaOut = 1
    lpaout(1) = cataOptpara(cataNbIn+ipara)
    call rsexch(' ', resultOut, option, numeStore, fieldOut, ierd)
    lchout(1) = fieldOut

    call jedema()
!
end subroutine
