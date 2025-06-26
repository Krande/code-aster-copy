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

subroutine liscom(nomo, codarr, lischa, l_need_modelz)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lisico.h"
#include "asterfort/lislch.h"
#include "asterfort/lislco.h"
#include "asterfort/lisnnb.h"
#include "asterfort/utmess.h"
    character(len=19) :: lischa
    character(len=1) :: codarr
    character(len=8) :: nomo
    aster_logical, intent(in), optional :: l_need_modelz
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (LISTE_CHARGES)
!
! VERIFICATION DE LA COHERENCE ET DE LA PRESENCE DES MODELES
!
! ----------------------------------------------------------------------
!
!
! IN  NOMO   : NOM DU MODELE
! IN  CODARR : TYPE DE MESSAGE INFO/ALARME/ERREUR SI PAS COMPATIBLE
! IN  LISCHA : SD LISTE DES CHARGES
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ichar, nbchar
    character(len=8) :: charge, modch1
    integer(kind=8) :: genrec
    aster_logical :: lveag, lveas, l_need_model
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOMBRE DE CHARGES
!
    call lisnnb(lischa, nbchar)
    if (nbchar .eq. 0) goto 999

    if (present(l_need_modelz)) then
        l_need_model = l_need_modelz
    else
        l_need_model = ASTER_FALSE
    end if
!
! --- BOUCLE SUR LES CHARGES
!
    do ichar = 1, nbchar
        call lislch(lischa, ichar, charge)
        call lislco(lischa, ichar, genrec)
        lveag = lisico('VECT_ASSE_GENE', genrec)
        lveas = lisico('VECT_ASSE', genrec)
        if (.not. lveag .and. .not. lveas) then
            if (nomo .ne. ' ') then
                call dismoi('NOM_MODELE', charge, 'CHARGE', repk=modch1)
                if (modch1 .ne. nomo) then
                    call utmess(codarr, 'CHARGES5_5', sk=charge)
                end if
            elseif (l_need_model) then
                call utmess(codarr, 'CHARGES5_6')
            end if
        end if
    end do
!
999 continue
!
    call jedema()
end subroutine
