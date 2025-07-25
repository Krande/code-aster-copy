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
subroutine nmsigi(ligrmo, compor, sigini)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/calcul.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    character(len=24) :: ligrmo, compor
    character(len=19) :: sigini
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (INITIALISATION)
!
! RECUPERATION D'UN CHARGEMENT MECANIQUE DE TYPE PRECONTRAINTE
!
! ----------------------------------------------------------------------
!
!
! IN  LIGRMO : NOM DU LIGREL DU MODELE
! IN  COMPOR : NOM DE LA CARTE DE COMPORTEMENT
! I/O SIGINI : NOM DU CHAMP DE CONTRAINTE INITIAL VIERGE.
!                    ON LUI ADDITIONNE LES EVENTUELS CHAMPS :
!                    LCHAR(I)//.CHME.SIGIN'
!
!
!
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=2)
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    integer(kind=8) :: ibid, iocc, iret, nbocc
    character(len=8) :: charge
    character(len=16) :: option
    character(len=19) :: sigcha
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call getfac('EXCIT', nbocc)
    if (nbocc .gt. 0) then
        option = 'ADD_SIGM'
        lchin(1) = sigini
        lpain(1) = 'PEPCON1'
        lpain(2) = 'PEPCON2'
        lchout(1) = '&&NMSIGI.PEPCON3'
        lpaout(1) = 'PEPCON3'
        do iocc = 1, nbocc
            call getvid('EXCIT', 'CHARGE', iocc=iocc, scal=charge, nbret=ibid)
            sigcha = charge//'.CHME.SIGIN'
            call exisd('CHAMP_GD', sigcha, iret)
            if (iret .gt. 0) then
                lchin(2) = sigcha
                call copisd('CHAM_ELEM_S', 'V', compor, lchout(1))
                call calcul('S', option, ligrmo, 2, lchin, &
                            lpain, 1, lchout, lpaout, 'V', &
                            'OUI')
                call detrsd('CHAM_ELEM_S', lchout(1))
                call copisd('CHAMP_GD', 'V', '&&NMSIGI.PEPCON3', sigini)
            end if
        end do
        call detrsd('CHAMP_GD', '&&NMSIGI.PEPCON3')
    end if
!
    call jedema()
!.
end subroutine
