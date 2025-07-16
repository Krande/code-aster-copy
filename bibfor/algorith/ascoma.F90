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
subroutine ascoma(hval_meelem, numeDof, listLoad, matrAsse)
!
    implicit none
!
#include "asterfort/asmatr.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/reajre.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: hval_meelem(*)
    character(len=24), intent(in) :: numeDof
    character(len=19), intent(in) :: listLoad, matrAsse
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! ASSEMBLAGE DE LA MATRICE DE RIGIDITE ASSOCIEE AUX CHARGEMENTS
! SUIVEURS
!
! --------------------------------------------------------------------------------------------------
!
! In  hval_meelem : hat-variable for elementary matrices
! IN  MEELEM : LISTE DES MATR_ELEM
! IN  NUMEDD : NOM DE LA NUMEROTATION MECANIQUE
! IN  LISCHA : SD L_CHARGE
! OUT MATASS : MATRICE GLOBALE ASSEMBLEE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbCoef, iret, iCoef
    character(len=24) :: coefJvName
    character(len=19) :: mesuiv
    character(len=24), pointer :: relr(:) => null()
    real(kind=8), pointer :: listCoef(:) => null(), coefMatr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Get elementary matrices for undead load
    call nmchex(hval_meelem, 'MEELEM', 'MESUIV', mesuiv)
    coefJvName = mesuiv(1:15)//'.COEF'

    call jeexin(coefJvName, iret)
    if (iret .ne. 0) then
        call jelira(coefJvName, 'LONUTI', nbCoef)
        call jeveuo(mesuiv(1:19)//'.RELR', 'L', vk24=relr)
        call jeveuo(coefJvName, 'L', vr=listCoef)
        call jedupo(mesuiv(1:19)//'.RERR', 'V', '&&ASCOMA           .RERR', ASTER_TRUE)
        call wkvect('&&ASCOMA.LISTE_COEF', 'V V R', 1, vr=coefMatr)
        do iCoef = 1, nbCoef
            call jedetr('&&ASCOMA           .RELR')
            call reajre('&&ASCOMA', relr(iCoef), 'V')
            coefMatr(1) = listCoef(iCoef)
            call asmatr(1, '&&ASCOMA           ', '&&ASCOMA.LISTE_COEF', numeDof, listLoad, &
                        'CUMU', 'V', 1, matrAsse)
        end do
        call jedetr('&&ASCOMA           .RELR')
        call jedetr('&&ASCOMA           .RERR')
        call jedetr('&&ASCOMA.LISTE_COEF')
    end if
!
    call jedema()
end subroutine
