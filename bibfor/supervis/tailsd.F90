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

subroutine tailsd(nom, nomsd, val, nbval)
    implicit none
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/rsorac.h"
    integer(kind=8) :: nbval, val(nbval)
    character(len=*) :: nom, nomsd
! ---------------------------------------------------------------
!  DETERMINE LE NOMBRE DE VALEURS DANS UNE STRUCTURE DE DONNES
! ---------------------------------------------------------------
! IN  NOM     K*  NOM IDENTIFIANT A LA FOIS LA SD ET LA METHODE
! IN  NOMSD   I   NOM DE LA SD A INTERROGER
! OUT VAL     I   VECTEUR D'ENTIERS DE TAILLE NBVAL CONTENANT LES
!                 DIFFERENTES TAILLES
! IN  NBVAL   I   TAILLE DU VECTEUR D'ENTIERS VAL
! ---------------------------------------------------------------
    character(len=19) :: sd19
    character(len=24) :: sd
    integer(kind=8) :: iret1, iret2, ibid
    complex(kind=8) :: cbid
!
!  DETERMINE LE NOMBRE MAXIMUM DE CHAMPS ET DE PARAMETRES
!  AINSI QUE LE NOMBRE EFFECTIF DE NUMEROS D'ORDRE POUR UNE SD RESULTAT
! ---------------------------------------------------------------------
    if (nom .eq. 'LIST_RESULTAT') then
        val(1) = 0
!       0 EST UTILISE POUR LES NUMEROS D'ORDRE (CF. RSACPA)
        val(2) = -1
        val(3) = 0
!
        sd19 = nomsd
        call jeexin(sd19//'.DESC', iret1)
        call jeexin(sd19//'.NOVA', iret2)
        if (iret1 .eq. 0 .or. iret2 .eq. 0) then
            goto 999
        end if
        call jelira(sd19//'.DESC', 'NOMMAX', val(1))
        call jelira(sd19//'.NOVA', 'NOMMAX', val(2))
        call rsorac(sd19, 'LONUTI', 0, 0.d0, ' ', &
                    cbid, 0.d0, ' ', val(3), 1, &
                    ibid)
!
!
!  DETERMINE LE NOMBRE D OBJETS SIMPLES DANS UNE COLECTION
! ---------------------------------------------------------------------
    else if (nom .eq. 'LIST_COLLECTION') then
        val(1) = 0
!
        sd = nomsd
        call jeexin(sd, iret1)
        if (iret1 .eq. 0) then
            goto 999
        end if
        call jelira(sd, 'NUTIOC', val(1))
        if (val(1) .eq. 0) call jelira(sd, 'NMAXOC', val(1))
    end if
999 continue
end subroutine
